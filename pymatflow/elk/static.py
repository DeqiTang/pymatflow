#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import copy
import shutil

from pymatflow.remote.server import server_handle
from pymatflow.octopus.octopus import octopus

"""
usage:
"""

class static_run(octopus):
    """
    """
    def __init__(self):
        super().__init__()

        #self.inp.set_runtype(runtype="static")
        self.magnetic_status = "spin-unpolarized" # "spin-polarized" "non-collinear"

    def scf(self, directory="tmp-octopus-static", runopt="gen", auto=0):
        """
        directory: a place for all the generated files

        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it

        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            #shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.inp[0].system.xyz.file, directory))

            #with open(os.path.join(directory, "POSCAR"), 'w') as fout:
            #    self.poscar.to_poscar(fout)

            # gen llhpc script
            self.gen_llhpc(directory=directory, scriptname="static-scf.slurm", cmd="$PMF_OCTOPUS")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="$PMF_OCTOPUS", scriptname="static-scf.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen local bash script
            self.gen_bash(directory=directory, mpi=self.run_params["mpi"], cmd="$PMF_OCTOPUS", scriptname="static-scf.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="$PMF_OCTOPUS", scriptname="static-scf.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"])

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("vasp")
            os.system("bash static-scf.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static-scf", server=self.run_params["server"])


    def nscf(self, directory="tmp-octopus-static", runopt="gen", auto=0):
        """
        directory: a place for all the generated files

        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("non-scf calculation:\n")
            print("  directory of previous scf calculattion not found!\n")
            sys.exit(1)

        if runopt == "gen" or runopt == "genrun":

            # gen llhpc script
            self.gen_llhpc(directory=directory, scriptname="static-nscf.slurm", cmd="$PMF_OCTOPUS")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="$PMF_OCTOPUS", scriptname="static-nscf.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # gen local bash script
            self.gen_bash(directory=directory, mpi=self.run_params["mpi"], cmd="$PMF_OCTOPUS", scriptname="static-nscf.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="$PMF_OCTOPUS", scriptname="static-nscf.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"])


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("vasp")
            os.system("bash static-nscf.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static-nscf", server=self.run_params["server"])

    def bands(self, directory="tmp-octopus-static", runopt="gen", auto=0):
        """
        directory: a place for all the generated files

        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it
        """
        self.set_params({
            "Calculation Modes/CalculationMode": "unocc",
            })
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("bands calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)

        if runopt == "gen" or runopt == "genrun":

            # gen llhpc script
            self.gen_llhpc(directory=directory, scriptname="static-bands.slurm", cmd="$PMF_OCTOPUS")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="$PMF_OCTOPUS", scriptname="static-bands.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # gen local bash script
            self.gen_bash(directory=directory, cmd="%s %PMF_OCTOPUS" % self.run_params["mpi"], scriptname="static-bands.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="$PMF_OCTOPUS", scriptname="static-bands.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"])

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("vasp")
            os.system("bash static-bands.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static-bands", server=self.run_params["server"])

    def converge_spacing(self, spacing_range=[], directory="tmp-octopus-encut", runopt="gen", auto=0, restart=0):
        """
        spacing_range: [[spacing_x, spacing_y, spacing_z], [spacing_x, spacing_y, spacing_z],...]
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.poscar.xyz.file, directory))

            os.chdir(directory)
            n_test = int((emax - emin) / step)
            for i in range(n_test + 1):
                encut = int(emin + i * step)
                os.mkdir("encut-%d" % encut)
                self.incar.params['ENCUT'] = encut
                with open(os.path.join("encut-%d" % encut, "INCAR"), 'w') as fout:
                    self.incar.to_incar(fout)
                with open(os.path.join("encut-%d" % encut, "POSCAR"), 'w') as fout:
                    self.poscar.to_poscar(fout)
                with open(os.path.join("encut-%d" % encut, "KPOINTS"), 'w') as fout:
                    self.kpoints.to_kpoints(fout)
                shutil.copyfile("POTCAR", os.path.join("encut-%d" % encut, "POTCAR"))


            # gen llhpc script
            with open("converge-encut.slurm", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                for i in range(n_test + 1):
                    encut = int(emin + i * step)
                    fout.write("cd ./encut-%d\n" % encut)
                    fout.write("yhrun $PMF_OCTOPUS\n")
                    fout.write("cd ../\n")
                    fout.write("\n")
            # gen pbs running script
            with open("converge-encut.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    encut = int(emin + i * step)
                    fout.write("cd ./encut-%d\n" % encut)
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi %s\n" % ("$PMF_OCTOPUS"))
                    fout.write("cd ../\n")
                    fout.write("\n")
            os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            for i in range(n_test + 1):
                encut = int(emin + i * step)
                os.chdir("encut-%d" % encut)
                os.system("%s vasp" % (self.run_params["mpi"]))
                os.chdir("../")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="converge-encut", server=self.run_params["server"])

    def converge_kpoints(self, kmin, kmax, step, directory="tmp-octopus-kpoints", runopt="gen", auto=0):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.poscar.xyz.file, directory))


            os.chdir(directory)
            n_test = int((kmax - kmin) / step)
            for i in range(n_test + 1):
                kpoints = int(kmin + i * step)
                os.mkdir("kpoints-%d" % kpoints)
                self.kpoints.set_kpoints(kpoints_mp=[kpoints, kpoints, kpoints, 0, 0, 0])
                with open(os.path.join("kpoints-%d" % kpoints, "INCAR"), 'w') as fout:
                    self.incar.to_incar(fout)
                with open(os.path.join("kpoints-%d" % kpoints, "POSCAR"), 'w') as fout:
                    self.poscar.to_poscar(fout)
                with open(os.path.join("kpoints-%d" % kpoints, "KPOINTS"), 'w') as fout:
                    self.kpoints.to_kpoints(fout)
                shutil.copyfile("POTCAR", os.path.join("kpoints-%d" % kpoints, "POTCAR"))


            # gen llhpc script
            with open("converge-kpoints.slurm", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                for i in range(n_test + 1):
                    kpoints = int(kmin + i * step)
                    fout.write("cd ./kpoints-%d\n" % kpoints)
                    fout.write("yhrun $PMF_OCTOPUS\n")
                    fout.write("cd ../\n")
                    fout.write("\n")
            # gen pbs running script
            with open("converge-kpoints.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    kpoints = int(kmin + i * step)
                    fout.write("cd ./kpoints-%d\n" % kpoints)
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi %s\n" % ("$PMF_OCTOPUS"))
                    fout.write("cd ../\n")
                    fout.write("\n")
            os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            for i in range(n_test + 1):
                kpoints = int(kmin + i * step)
                os.chdir("kpoints-%d" % kpoints)
                os.system("%s vasp" % (self.run_params["mpi"]))
                os.chdir("../")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="converge-kpoints", server=self.run_params["server"])

    def converge_sigma(self, sigma_min, sigma_max, step, directory="tmp-octopus-sigma", runopt="gen", auto=0):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.poscar.xyz.file, directory))


            os.chdir(directory)
            n_test = int((sigma_max - sigma_min) / step)
            for i in range(n_test + 1):
                sigma = sigma_min + i * step
                os.mkdir("sigma-%.6f" % sigma)
                self.incar.electrons.params['SIGMA'] = sigma
                with open(os.path.join("sigma-%.6f" % sigma, "INCAR"), "w") as fout:
                    self.incar.to_incar(fout)
                with open(os.path.join("sigma-%.6f" % sigma, "POSCAR"), 'w') as fout:
                    self.poscar.to_poscar(fout)
                with open(os.path.join("sigma-%.6f" % sigma, "KPOINTS"), 'w') as fout:
                    self.kpoints.to_kpoints(fout)
                shutil.copyfile("POTCAR", os.path.join("sigma-%.6f" % sigma, "POTCAR"))


            # gen llhpc running script
            with open("converge-sigma.slurm", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                for i in range(n_test + 1):
                    sigma = sigma_min + i * step
                    fout.write("cd ./sigma-%.6f\n" % sigma)
                    fout.write("yhrun %PMF_OCTOPUS\n")
                    fout.write("cd ../\n")
                    fout.write("\n")
            # gen pbs running script
            with open("converge-sigma.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % (self.run_params["jobname"]))
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    sigma = sigma_min + i * step
                    fout.write("cd ./sigma-%.6f\n" % sigma)
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi %s\n" % ("$PMF_OCTOPUS"))
                    fout.write("cd ../\n")
                    fout.write("\n")
            os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            for i in range(n_test + 1):
                sigma = sigma_min + i * step
                os.chdir("sigma-%.6f" % sigma)
                os.system("%s vasp" % (self.run_params["mpi"]))
                os.chdir("../")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="converge-sigma", server=self.run_params["server"])

    #

    def set_scf(self, params):
        pass

    def band(self, directory="tmp-octopus-static", runopt="gen", auto=0, kpath=None):
        """
        directory: a place for all the generated files
        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it
        Note: scf nscf(bands) in a single run
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            #shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.inp[0].system.xyz.file, directory))

            # scf
            self.inp[0].set_params({
                "Calculation Modes/CalculationMode": "gs",
            })


            # nscf: bands
            self.inp[1].set_params({
                "Calculation Modes/CalculationMode": "unocc",
                })            

            # gen llhpc script
            with open(os.path.join(directory, "static.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                
                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("yhrun $PMF_OCTOPUS > output.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")                
                fout.write("yhrun $PMF_OCTOPUS > output.nscf\n")


            # gen pbs script
            with open(os.path.join(directory, "static.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_OCTOPUS > output.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_OCTOPUS > output.nscf\n")                
                

            # gen local bash script
            with open(os.path.join(directory, "static.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("%s $PMF_OCTOPUS > output.scf\n" % self.run_params["mpi"])
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("%s $PMF_OCTOPUS > output.nscf\n" % self.run_params["mpi"])
                

            # gen lsf_sz script
            with open(os.path.join(directory, "static.lsf_sz"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("APP_NAME=%s\n" % self.run_params["queue"])
                fout.write("NP=%d\n" % (self.run_params["nodes"] * self.run_params["ppn"]))
                fout.write("NP_PER_NODE=%d\n" % self.run_params["ppn"])
                fout.write("RUN=\"RAW\"\n")
                fout.write("CURDIR=$PWD\n")
                fout.write("#VASP=/home-yg/Soft/Vasp5.4/vasp_std\n")
                fout.write("source /home-yg/env/intel-12.1.sh\n")
                fout.write("source /home-yg/env/openmpi-1.6.5-intel.sh\n")
                fout.write("cd $CURDIR\n")
                fout.write("# starting creating ./nodelist\n")
                fout.write("rm -rf $CURDIR/nodelist >& /dev/null\n")
                fout.write("for i in `echo $LSB_HOSTS`\n")
                fout.write("do\n")
                fout.write("  echo \"$i\" >> $CURDIR/nodelist \n")
                fout.write("done\n")
                fout.write("ndoelist=$(cat $CURDIR/nodelist | uniq | awk \'{print $1}\' | tr \'\n\' \',\')\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_OCTOPUS > output.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_OCTOPUS > output.nscf\n")


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash static.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static", server=self.run_params["server"])


    def dos(self, directory="tmp-octopus-static", runopt="gen", auto=0):
        """
        directory: a place for all the generated files

        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it
        Note: scf nscf(dos) in a single run
        """
        """
        directory: a place for all the generated files
        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it
        Note: scf nscf(bands) in a single run
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            #shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.inp[0].system.xyz.file, directory))

            # scf
            self.inp[0].set_params({
                "Calculation Modes/CalculationMode": "gs",
            })


            # nscf: bands
            self.inp[1].set_params({
                "Calculation Modes/CalculationMode": "unocc",
                })            
            

            # gen llhpc script
            with open(os.path.join(directory, "static.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                
                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("yhrun $PMF_OCTOPUS > output.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")                
                fout.write("yhrun $PMF_OCTOPUS \n")


            # gen pbs script
            with open(os.path.join(directory, "static.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_OCTOPUS > output.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_OCTOPUS > output.nscf\n")                
                

            # gen local bash script
            with open(os.path.join(directory, "static.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("%s $PMF_OCTOPUS > output.scf\n" % self.run_params["mpi"])
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("%s $PMF_OCTOPUS > output.nscf\n" % self.run_params["mpi"])
                

            # gen lsf_sz script
            with open(os.path.join(directory, "static.lsf_sz"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("APP_NAME=%s\n" % self.run_params["queue"])
                fout.write("NP=%d\n" % (self.run_params["nodes"] * self.run_params["ppn"]))
                fout.write("NP_PER_NODE=%d\n" % self.run_params["ppn"])
                fout.write("RUN=\"RAW\"\n")
                fout.write("CURDIR=$PWD\n")
                fout.write("#VASP=/home-yg/Soft/Vasp5.4/vasp_std\n")
                fout.write("source /home-yg/env/intel-12.1.sh\n")
                fout.write("source /home-yg/env/openmpi-1.6.5-intel.sh\n")
                fout.write("cd $CURDIR\n")
                fout.write("# starting creating ./nodelist\n")
                fout.write("rm -rf $CURDIR/nodelist >& /dev/null\n")
                fout.write("for i in `echo $LSB_HOSTS`\n")
                fout.write("do\n")
                fout.write("  echo \"$i\" >> $CURDIR/nodelist \n")
                fout.write("done\n")
                fout.write("ndoelist=$(cat $CURDIR/nodelist | uniq | awk \'{print $1}\' | tr \'\n\' \',\')\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_OCTOPUS > output.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_OCTOPUS > output.nscf\n")


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash static.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static", server=self.run_params["server"])


    def optics(self, directory="tmp-octopus-static", runopt="gen", auto=0):
        """
        directory: a place for all the generated files

        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it
        Note: scf nscf(pdos, bands) in a single run
        """
        """
        directory: a place for all the generated files
        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it
        Note: scf nscf(bands) in a single run
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            #shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.inp[0].system.xyz.file, directory))

            # scf
            self.inp[0].set_params({
                "Calculation Modes/CalculationMode": "gs",
            })


            # nscf: bands
            self.inp[1].set_params({
                "Calculation Modes/CalculationMode": "unocc",
                })            

            # gen llhpc script
            with open(os.path.join(directory, "static.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                
                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("yhrun $PMF_OCTOPUS > output.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")                
                fout.write("yhrun $PMF_OCTOPUS > output.nscf\n")


            # gen pbs script
            with open(os.path.join(directory, "static.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_OCTOPUS > output.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_OCTOPUS > output.nscf\n")                
                

            # gen local bash script
            with open(os.path.join(directory, "static.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("%s $PMF_OCTOPUS > output.scf\n" % self.run_params["mpi"])
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("%s $PMF_OCTOPUS > output.nscf\n" % self.run_params["mpi"])
                

            # gen lsf_sz script
            with open(os.path.join(directory, "static.lsf_sz"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("APP_NAME=%s\n" % self.run_params["queue"])
                fout.write("NP=%d\n" % (self.run_params["nodes"] * self.run_params["ppn"]))
                fout.write("NP_PER_NODE=%d\n" % self.run_params["ppn"])
                fout.write("RUN=\"RAW\"\n")
                fout.write("CURDIR=$PWD\n")
                fout.write("#VASP=/home-yg/Soft/Vasp5.4/vasp_std\n")
                fout.write("source /home-yg/env/intel-12.1.sh\n")
                fout.write("source /home-yg/env/openmpi-1.6.5-intel.sh\n")
                fout.write("cd $CURDIR\n")
                fout.write("# starting creating ./nodelist\n")
                fout.write("rm -rf $CURDIR/nodelist >& /dev/null\n")
                fout.write("for i in `echo $LSB_HOSTS`\n")
                fout.write("do\n")
                fout.write("  echo \"$i\" >> $CURDIR/nodelist \n")
                fout.write("done\n")
                fout.write("ndoelist=$(cat $CURDIR/nodelist | uniq | awk \'{print $1}\' | tr \'\n\' \',\')\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_OCTOPUS > outpu.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_OCTOPUS > output.nscf\n")


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash static.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static", server=self.run_params["server"])


    def bse(self, directory="tmp-octopus-static", runopt="gen", auto=0, bse_level=0, algo_gw="EVGW"):
        """
        directory: a place for all the generated files

        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it
        bse_level:    
            0 -> bse on standard DFT
            1 -> bse on hybrid functional 
            2 -> bse on GW            
        Note:
        Reference:
            https://www.vasp.at/wiki/index.php/Practical_guide_to_GW_calculations
            https://www.vasp.at/wiki/index.php/BSE_calculations
            https://www.vasp.at/wiki/index.php/Dielectric_properties_of_Si_using_BSE
            https://www.vasp.at/wiki/index.php/Plotting_the_BSE_fatband_structure_of_Si
            https://www.vasp.at/wiki/index.php/Bandgap_of_Si_in_GW
            https://www.vasp.at/wiki/index.php/Bandstructure_of_Si_in_GW_(VASP2WANNIER90)
            https://www.vasp.at/wiki/index.php/Bandstructure_of_SrVO3_in_GW
            https://www.vasp.at/wiki/index.php/Improving_the_dielectric_function
        """
        """
        directory: a place for all the generated files
        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it
        Note: scf nscf(bands) in a single run
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            #shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.inp[0].system.xyz.file, directory))

            # scf
            self.inp[0].set_params({
                "Calculation Modes/CalculationMode": "gs",
            })


            # nscf: bands
            self.inp[1].set_params({
                "Calculation Modes/CalculationMode": "unocc",
                })            

            # gen llhpc script
            with open(os.path.join(directory, "static.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                
                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("yhrun $PMF_OCTOPUS > output.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")                
                fout.write("yhrun $PMF_OCTOPUS > output.nscf\n")


            # gen pbs script
            with open(os.path.join(directory, "static.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_OCTOPUS > output.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_OCTOPUS > output.nscf\n")                
                

            # gen local bash script
            with open(os.path.join(directory, "static.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("%s $PMF_OCTOPUS > output.scf\n" % self.run_params["mpi"])
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("%s $PMF_OCTOPUS > output.nscf\n" % self.run_params["mpi"])
                

            # gen lsf_sz script
            with open(os.path.join(directory, "static.lsf_sz"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("APP_NAME=%s\n" % self.run_params["queue"])
                fout.write("NP=%d\n" % (self.run_params["nodes"] * self.run_params["ppn"]))
                fout.write("NP_PER_NODE=%d\n" % self.run_params["ppn"])
                fout.write("RUN=\"RAW\"\n")
                fout.write("CURDIR=$PWD\n")
                fout.write("#VASP=/home-yg/Soft/Vasp5.4/vasp_std\n")
                fout.write("source /home-yg/env/intel-12.1.sh\n")
                fout.write("source /home-yg/env/openmpi-1.6.5-intel.sh\n")
                fout.write("cd $CURDIR\n")
                fout.write("# starting creating ./nodelist\n")
                fout.write("rm -rf $CURDIR/nodelist >& /dev/null\n")
                fout.write("for i in `echo $LSB_HOSTS`\n")
                fout.write("do\n")
                fout.write("  echo \"$i\" >> $CURDIR/nodelist \n")
                fout.write("done\n")
                fout.write("ndoelist=$(cat $CURDIR/nodelist | uniq | awk \'{print $1}\' | tr \'\n\' \',\')\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_OCTOPUS > output.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_OCTOPUS > output.nscf\n")


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash static.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static", server=self.run_params["server"])


    def stm(self, directory="tmp-octopus-static", runopt="gen", auto=0, hse_in_scf=True):
        """
        directory: a place for all the generated files

        hse_in_scf:
            if HSE is used, choose whether to use HSE in both scf and nscf or only in nscf.
            if HSE is not used, hse_in_scf will do nothing
        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it
        Note: scf stm in a single run
        """
        """
        directory: a place for all the generated files
        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it
        Note: scf nscf(bands) in a single run
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            #shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.inp[0].system.xyz.file, directory))

            # scf
            self.inp[0].set_params({
                "Calculation Modes/CalculationMode": "gs",
            })


            # nscf: bands
            self.inp[1].set_params({
                "Calculation Modes/CalculationMode": "unocc",
                })            

            # gen llhpc script
            with open(os.path.join(directory, "static.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                
                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("yhrun $PMF_OCTOPUS > output.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")                
                fout.write("yhrun $PMF_OCTOPUS > output.nscf\n")


            # gen pbs script
            with open(os.path.join(directory, "static.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_OCTOPUS > output.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_OCTOPUS > output.nscf\n")                
                

            # gen local bash script
            with open(os.path.join(directory, "static.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("%s $PMF_OCTOPUS > output.scf\n" % self.run_params["mpi"])
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("%s $PMF_OCTOPUS > output.nscf\n" % self.run_params["mpi"])
                

            # gen lsf_sz script
            with open(os.path.join(directory, "static.lsf_sz"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("APP_NAME=%s\n" % self.run_params["queue"])
                fout.write("NP=%d\n" % (self.run_params["nodes"] * self.run_params["ppn"]))
                fout.write("NP_PER_NODE=%d\n" % self.run_params["ppn"])
                fout.write("RUN=\"RAW\"\n")
                fout.write("CURDIR=$PWD\n")
                fout.write("#VASP=/home-yg/Soft/Vasp5.4/vasp_std\n")
                fout.write("source /home-yg/env/intel-12.1.sh\n")
                fout.write("source /home-yg/env/openmpi-1.6.5-intel.sh\n")
                fout.write("cd $CURDIR\n")
                fout.write("# starting creating ./nodelist\n")
                fout.write("rm -rf $CURDIR/nodelist >& /dev/null\n")
                fout.write("for i in `echo $LSB_HOSTS`\n")
                fout.write("do\n")
                fout.write("  echo \"$i\" >> $CURDIR/nodelist \n")
                fout.write("done\n")
                fout.write("ndoelist=$(cat $CURDIR/nodelist | uniq | awk \'{print $1}\' | tr \'\n\' \',\')\n")

                fout.write("#scf\n")
                fout.write("cat >inp<<EOF\n")
                fout.write(self.inp[0].to_string())
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_OCTOPUS > output.scf\n")
                
                fout.write("# nscf\n")
                fout.write("cat > inp<<EOF\n")
                fout.write(self.inp[1].to_string())
                fout.write("EOF\n")        
                fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_OCTOPUS > output.nscf\n")


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash static.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static", server=self.run_params["server"])