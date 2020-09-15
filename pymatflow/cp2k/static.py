"""
static calculation
"""
import os
import sys
import shutil
import numpy as np

from pymatflow.remote.server import server_handle
from pymatflow.cp2k.cp2k import cp2k

"""
"""

class static_run(cp2k):
    """
    Usage:
        a = static_run()
        a.get_xyz(xxx)
        a.set_params(xxx)
        a.set_printout(xxx)
        a.scf(xxx)
    Note:
        static_run is the class as an agent for static type calculation, including
        scf, CUTOFF converge test, REL_CUTOFF converge test. and it can control the
        calculation of properties like electronic band structure, projected density
        of states(pdos), electron density, elf, charge difference, etc.
    """
    def __init__(self):
        """
        TODO:
            include implement MP2 calculation through CP2K/ATOM
        """
        super().__init__()
        #self.glob = cp2k_glob()
        #self.force_eval = cp2k_force_eval()

        self.glob.basic_setting(run_type="ENERGY_FORCE")
        self.force_eval.basic_setting()


    def scf(self, directory="tmp-cp2k-static", inpname="static-scf.inp", output="static-scf.out", runopt="gen", auto=0):
        """
        :param directory:
            directory is and path where the calculation will happen.
        :param inpname:
            input filename for the cp2k
        :param output:
            output filename for the cp2k
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)

            # gen server job comit file
            self.gen_llhpc(directory=directory, inpname=inpname, output=output, cmd="$PMF_CP2K")
            # gen pbs server job comit file
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_CP2K", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])

        if runopt == "run" or runopt == "genrun":
           os.chdir(directory)
           os.system("%s $PMF_CP2K -in %s | tee %s" % (self.run_params["mpi"], inpname, output))
           os.chdir("../")

        server_handle(auto=auto, directory=directory, jobfilebase="static-scf", server=self.run_params["server"])

    def converge_cutoff(self, emin, emax, step, directory="tmp-cp2k-cutoff", runopt="gen", auto=0):
        """
        Note:
            this function is used to do the converge test for CUTOFF.
            advices on CUTOFF converging:
            we can first check the basis set file, and find the largest
            value among all elements used. and we times it by 4. then
            we can set the converge range convering that value, to find
            the converged CUTOFF
        :param emin:
            the minimum cutoff of the test range
        :param emax:
            the maximum cutoff of the test range
        :param step:
            the step for the converge test
        :param rel_cutoff:
            for the test of cutoff, rel_cutoff is set to an fixed value
        :param directory:
            where the converge test happens
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

            n_test = int((emax - emin) / step)
            for i in range(n_test + 1):
                cutoff = int(emin + i * step)
                inpname = "cutoff-%d.inp" % cutoff
                self.force_eval.dft.mgrid.params["CUTOFF"] = cutoff
                #self.force_eval.dft.mgrid.params["REL_CUTOFF"] = rel_cutoff

                with open(os.path.join(directory, inpname), 'w') as fout:
                    self.glob.to_input(fout)
                    self.force_eval.to_input(fout)

            # gen yhbatch running script
            with open(os.path.join(directory, "converge-cutoff.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                for i in range(n_test + 1):
                    cutoff = int(emin + i * step)
                    inpname = "cutoff-%d.inp" % cutoff
                    out_f_name = "cutoff-%d.out" % cutoff
                    fout.write("yhrun $PMF_CP2K -in %s > %s\n" % (inpname, out_f_name))

            # gen pbs running script
            with open(os.path.join(directory, "converge-cutoff.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    cutoff = int(emin + i * step)
                    inpname = "cutoff-%d.inp" % cutoff
                    out_f_name = "cutoff-%d.out" % cutoff
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_CP2K -in %s > %s\n" % (inpname, out_f_name))

        # run
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for i in range(n_test + 1):
                cutoff = int(emin + i * step)
                inpname = "cutoff-%d.inp" % cutoff
                output = "cutoff-%d.out" % cutoff
                os.system("%s $PMF_CP2K -in %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="converge-cutoff", server=self.run_params["server"])

    def converge_rel_cutoff(self, emin, emax, step, directory="tmp-cp2k-rel-cutoff", runopt="gen", auto=0):
        """
        Note:
            this function is used to do the converge test of REL_CUTOFF.
        :param emin:
            the minimum rel_cutoff of the test range
        :param emax:
            the maximum rel_cutoff of the test range
        :param step:
            the step for the converge test
        :param cutoff:
            for the test of rel_cutoff, cutoff is set to an fixed value
        :param directory:
            where the converge test happens
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

            n_test = int((emax - emin) / step)
            for i in range(n_test + 1):
                rel_cutoff = int(emin + i * step)
                inpname = "rel-cutoff-%d.inp" % rel_cutoff
                #self.force_eval.dft.mgrid.params["CUTOFF"] = cutoff
                self.force_eval.dft.mgrid.params["REL_CUTOFF"] = rel_cutoff

                with open(os.path.join(directory, inpname), 'w') as fout:
                    self.glob.to_input(fout)
                    self.force_eval.to_input(fout)

            # gen yhbatch running script
            with open(os.path.join(directory, "converge-rel-cutoff.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                for i in range(n_test + 1):
                    rel_cutoff = int(emin + i * step)
                    inpname = "rel-cutoff-%d.inp" % rel_cutoff
                    out_f_name = "rel-cutoff-%d.out" % rel_cutoff
                    fout.write("yhrun $PMF_CP2K -in %s > %s\n" % (inpname, out_f_name))

            # gen pbs running script
            with open(os.path.join(directory, "converge-rel-cutoff.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    rel_cutoff = int(emin + i * step)
                    inpname = "rel-cutoff-%d.inp" % rel_cutoff
                    out_f_name = "rel-cutoff-%d.out" % rel_cutoff
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_CP2K -in %s > %s\n" % (inpname, out_f_name))

        # run
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for i in range(n_test + 1):
                rel_cutoff = int(emin + i * step)
                inpname = "rel-cutoff-%d.inp" % rel_cutoff
                output = "rel-cutoff-%d.out" % rel_cutoff
                os.system("%s $PMF_CP2K -in %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="converge-rel-cutoff", server=self.run_params["server"])

    def converge_kpoints_auto(self, kmin, kmax, step, directory="tmp-cp2k-kpoints-auto", runopt="gen", auto=0):
        """
        Note:
            this function is used to do the converge test for KPONTS-AUTO.
            in this test mode, we input kmin, kmax, and step. the code will
            generate continuously the kpoint to test, like when kmin=1, kmax=3
            step=1. the tested kpoints would be: 1x1x1, 2x2x2, 3x3x3.
        :param kmin:
            the minimum kpoint of the test range
        :param kmax:
            the maximum kpoint of the test range
        :param step:
            the step for the converge test
        :param directory:
            where the converge test happens
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

            n_test = int((kmax - kmin) / step)
            for i in range(n_test + 1):
                kpoint = int(kmin + i * step)
                inpname = "kpoints-%d.inp" % kpoint

                self.force_eval.set_params({
                    "DFT-KPOINTS-SCHEME": "MONKHORST-PACK %d %d %d" % (kpoint, kpoint, kpoint)
                    })
                with open(os.path.join(directory, inpname), 'w') as fout:
                    self.glob.to_input(fout)
                    self.force_eval.to_input(fout)

            # gen yhbatch running script
            with open(os.path.join(directory, "converge-kpoints.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                for i in range(n_test + 1):
                    kpoint = int(kmin + i * step)
                    inpname = "kpoints-%d.inp" % kpoint
                    out_f_name = "kpoints-%d.out" % kpoint
                    fout.write("yhrun $PMF_CP2K -in %s > %s\n" % (inpname, out_f_name))

            # gen pbs running script
            with open(os.path.join(directory, "converge-kpoints.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    kpoint = int(kmin + i * step)
                    inpname = "kpoints-%d.inp" % kpoint
                    out_f_name = "kpoints-%d.out" % kpoint
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_CP2K -in %s > %s\n" % (inpname, out_f_name))

        # run
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for i in range(n_test + 1):
                kpoint = int(kmin + i * step)
                inpname = "kpoints-%d.inp" % kpoint
                output = "kpoints-%d.out" % kpoint
                os.system("%s $PMF_CP2K -in %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="converge-kpoints", server=self.run_params["server"])

    def converge_kpoints_manual(self, directory="tmp-cp2k-kpoints-manual", runopt="gen", auto=0,
            kpoints_list=[[1, 1, 1], [2, 2, 2], [3, 3, 3]]):
        """
        Note:
            this function is used to do the converge test for KPOINTS-MANUAL.
            in this mode, we have to specify clearly every k point to test, like
            kpoints_list=[[1, 1, 1], [1, 2, 1], [2, 2, 2]]
        :param: directory:
            where the converge test happens
        :param force_eval:
            allowing control of FORCE_EVAL/... parameters by user
        :param kpoints_list:
            kpoints test range
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

            for i in range(len(kpoints_list)):
                inpname = "kpoints-%dx%dx%d.inp" % (kpoints_list[i][0], kpoints_list[i][1], kpoints_list[i][2])
                self.force_eval.set_params({
                    "DFT-KPOINTS-SCHEME": "MONKHORST-PACK %d %d %d" % (kpoints_list[i][0], kpoints_list[i][1], kpoints_list[i][2])
                    })
                with open(os.path.join(directory, inpname), 'w') as fout:
                    self.glob.to_input(fout)
                    self.force_eval.to_input(fout)

            # gen yhbatch running script
            with open(os.path.join(directory, "converge-kpoints.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                for i in range(len(kpoints_list)):
                    inpname = "kpoints-%dx%dx%d.inp" % (kpoints_list[i][0], kpoints_list[i][1], kpoints_list[i][2])
                    out_f_name = "kpoints-%dx%dx%d.out" % (kpoints_list[i][0], kpoints_list[i][1], kpoints_list[i][2])
                    fout.write("yhrun $PMF_CP2K -in %s > %s\n" % (inpname, out_f_name))

            # gen obs running script
            with open(os.path.join(directory, "converge-kpoints.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(len(kpoints_list)):
                    inpname = "kpoints-%dx%dx%d.inp" % (kpoints_list[i][0], kpoints_list[i][1], kpoints_list[i][2])
                    out_f_name = "kpoints-%dx%dx%d.out" % (kpoints_list[i][0], kpoints_list[i][1], kpoints_list[i][2])
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_CP2K -in %s > %s\n" % (inpname, out_f_name))

        # run
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for i in range(len(kpoints_list)):
                inpname = "kpoints-%dx%dx%d.inp" % (kpoints_list[i][0], kpoints_list[i][1], kpoints_list[i][2])
                output = "kpoints-%dx%dx%d.out" % (kpoints_list[i][0], kpoints_list[i][1], kpoints_list[i][2])
                os.system("%s $PMF_CP2K -in %s > %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="converge-kpoints", server=self.run_params["server"])
