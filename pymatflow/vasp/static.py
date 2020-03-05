#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil

from pymatflow.remote.server import server_handle
from pymatflow.vasp.vasp import vasp

"""
usage:
"""

class static_run(vasp):
    """
    """
    def __init__(self):
        super().__init__()

        self.incar.set_runtype(runtype="static")


    def scf(self, directory="tmp-vasp-static", runopt="gen", auto=0):
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
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.poscar.xyz.file, directory))

            with open(os.path.join(directory, "POSCAR"), 'w') as fout:
                self.poscar.to_poscar(fout)

            # gen yhbatch script
            self.gen_yh(directory=directory, scriptname="static-scf.sub", cmd="vasp_std")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="vasp_std", scriptname="static-scf.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # gen local bash script
            self.gen_bash(directory=directory, mpi=self.run_params["mpi"], cmd="vasp_std", scriptname="static-scf.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="vasp_std", scriptname="static-scf.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"])

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("vasp")
            os.system("bash static-scf.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static-scf", server=self.run_params["server"])


    def nscf(self, directory="tmp-vasp-static", runopt="gen", auto=0):
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

            # gen yhbatch script
            self.gen_yh(directory=directory, scriptname="static-nscf.sub", cmd="vasp_std")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="vasp_std", scriptname="static-nscf.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # gen local bash script
            self.gen_bash(directory=directory, mpi=self.run_params["mpi"], cmd="vasp_std", scriptname="static-nscf.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="vasp_std", scriptname="static-nscf.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"])


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("vasp")
            os.system("bash static-nscf.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static-nscf", server=self.run_params["server"])

    def bands(self, directory="tmp-vasp-static", runopt="gen", auto=0):
        """
        directory: a place for all the generated files

        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it
        """
        self.set_params({
            "ICHARG": 11,
            "LORBIT": 11,
            })
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("bands calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)

        if runopt == "gen" or runopt == "genrun":

            # gen yhbatch script
            self.gen_yh(directory=directory, scriptname="static-bands.sub", cmd="vasp_std")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="vasp_std", scriptname="static-bands.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # gen local bash script
            self.gen_bash(directory=directory, cmd="%s vasp_std" % self.run_params["mpi"], scriptname="static-bands.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="vasp_std", scriptname="static-bands.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"])

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("vasp")
            os.system("bash static-bands.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static-bands", server=self.run_params["server"])

    def converge_encut(self, emin, emax, step, directory="tmp-vasp-encut", runopt="gen", auto=0, restart=0):
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


            # gen yhbatch running script
            with open("converge-encut.sub", 'w') as fout:
                fout.write("#!/bin/bash\n")
                for i in range(n_test + 1):
                    encut = int(emin + i * step)
                    fout.write("cd ./encut-%d\n" % encut)
                    fout.write("yhrun -N 1 -n 24 vasp\n")
                    fout.write("cd ../\n")
                    fout.write("\n")
            # gen pbs running script
            with open("converge-encut.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    encut = int(emin + i * step)
                    fout.write("cd ./encut-%d\n" % encut)
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi %s\n" % ("vasp_std"))
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

    def converge_kpoints(self, kmin, kmax, step, directory="tmp-vasp-kpoints", runopt="gen", auto=0):
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


            # gen yhbatch running script
            with open("converge-kpoints.sub", 'w') as fout:
                fout.write("#!/bin/bash\n")
                for i in range(n_test + 1):
                    kpoints = int(kmin + i * step)
                    fout.write("cd ./kpoints-%d\n" % kpoints)
                    fout.write("yhrun -N 1 -n 24 vasp\n")
                    fout.write("cd ../\n")
                    fout.write("\n")
            # gen pbs running script
            with open("converge-kpoints.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    kpoints = int(kmin + i * step)
                    fout.write("cd ./kpoints-%d\n" % kpoints)
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi %s\n" % ("vasp_std"))
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

    def converge_sigma(self, sigma_min, sigma_max, step, directory="tmp-vasp-sigma", runopt="gen", auto=0):
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


            # gen yhbatch running script
            with open("converge-sigma.sub", 'w') as fout:
                fout.write("#!/bin/bash\n")
                for i in range(n_test + 1):
                    sigma = sigma_min + i * step
                    fout.write("cd ./sigma-%.6f\n" % sigma)
                    fout.write("yhrun -N 1 -n 24 vasp\n")
                    fout.write("cd ../\n")
                    fout.write("\n")
            # gen pbs running script
            with open("converge-sigma.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % (self.run_params["jobname"]))
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    sigma = sigma_min + i * step
                    fout.write("cd ./sigma-%.6f\n" % sigma)
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi %s\n" % ("vasp_std"))
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

    def run(self, directory="tmp-vasp-static", runopt="gen", auto=0,
        kpoints_mp_scf=[1, 1, 1, 0, 0, 0], kpoints_mp_nscf=[3, 3, 3, 0, 0, 0], kpath=None, kpath_intersections=15):
        """
        directory: a place for all the generated files

        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it
        Note: scf nscf(pdos) bands in a single run
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.poscar.xyz.file, directory))

            # scf
            self.set_kpoints(option="automatic", kpoints_mp=kpoints_mp_scf)
            incar_scf = self.incar.to_string()
            kpoints_scf = self.kpoints.to_string()

            # nscf
            self.set_kpoints(option="automatic", kpoints_mp=kpoints_mp_nscf)
            incar_nscf = self.incar.to_string()
            kpoints_nscf = self.kpoints.to_string()

            # band structure
            self.set_params({
                "ICHARG": 11,
                "LORBIT": 11,
                })
            self.set_kpoints(option="bands", kpath=kpath, kpath_intersections=kpath_intersections)
            incar_bands = self.incar.to_string()
            kpoints_bands = self.kpoints.to_string()



            # gen yhbatch script
            with open(os.path.join(directory, "static.sub"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("cat >POSCAR<<EOF\n")
                self.poscar.to_poscar(fout)
                fout.write("EOF\n")
                fout.write("# scf\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_scf)
                fout.write("EOF\n")
                fout.write("cat > KPOINTS<<EOF\n")
                #self.kpoints.to_kpoints(fout)
                fout.write(kpoints_scf)
                fout.write("EOF\n")
                fout.write("yhrun -N 1 -n 24 $PMF_VASP_STD \n")

                fout.write("# nscf\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_nscf)
                fout.write("EOF\n")
                fout.write("cat >KPOINTS<<EOF\n")
                #self.kpoints.to_kpoints(fout)
                fout.write(kpoints_nscf)
                fout.write("EOF\n")
                fout.write("yhrun -N 1 -n 24 $PMF_VASP_STD \n")


                fout.write("# band structure\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_bands)
                fout.write("EOF\n")
                fout.write("cat >KPOINTS<<EOF\n")
                #self.kpoints.to_kpoints(fout)
                fout.write(kpoints_bands)
                fout.write("EOF\n")
                fout.write("yhrun -N 1 -n 24 $PMF_VASP_STD \n")

            # gen pbs script
            with open(os.path.join(directory, "static.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                fout.write("cat >POSCAR<<EOF\n")
                self.poscar.to_poscar(fout)
                fout.write("EOF\n")
                fout.write("# scf\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_scf)
                fout.write("EOF\n")
                fout.write("cat > KPOINTS<<EOF\n")
                #self.kpoints.to_kpoints(fout)
                fout.write(kpoints_scf)
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi $PMF_VASP_STD \n")

                fout.write("# nscf\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_nscf)
                fout.write("EOF\n")
                fout.write("cat >KPOINTS<<EOF\n")
                #self.kpoints.to_kpoints(fout)
                fout.write(kpoints_nscf)
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi $PMF_VASP_STD \n")


                fout.write("# band structure\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_bands)
                fout.write("EOF\n")
                fout.write("cat >KPOINTS<<EOF\n")
                #self.kpoints.to_kpoints(fout)
                fout.write(kpoints_bands)
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi $PMF_VASP_STD \n")

            # gen local bash script
            with open(os.path.join(directory, "static.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("cat >POSCAR<<EOF\n")
                self.poscar.to_poscar(fout)
                fout.write("EOF\n")
                fout.write("# scf\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_scf)
                fout.write("EOF\n")
                fout.write("cat > KPOINTS<<EOF\n")
                #self.kpoints.to_kpoints(fout)
                fout.write(kpoints_scf)
                fout.write("EOF\n")
                fout.write("%s $PMF_VASP_STD \n" % self.run_params["mpi"])

                fout.write("# nscf\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_nscf)
                fout.write("EOF\n")
                fout.write("cat >KPOINTS<<EOF\n")
                #self.kpoints.to_kpoints(fout)
                fout.write(kpoints_nscf)
                fout.write("EOF\n")
                fout.write("%s $PMF_VASP_STD \n" % self.run_params["mpi"])


                fout.write("# band structure\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_bands)
                fout.write("EOF\n")
                fout.write("cat >KPOINTS<<EOF\n")
                #self.kpoints.to_kpoints(fout)
                fout.write(kpoints_bands)
                fout.write("EOF\n")
                fout.write("%s $PMF_VASP_STD \n" % self.run_params["mpi"])

            # gen lsf_sz script
            with open(os.path.join(directory, "static.lsf_sz"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("APP_NAME=intelY_mid\n")
                fout.write("NP=%d\n" % (self.run_params["nodes"] * self.run_params["ppn"]))
                fout.write("NP_PER_NODE=%d\n" % self.run_params["ppn"])
                fout.write("RUN=\"RAW\"\n")
                fout.write("CURDIR=$PWD\n")
                fout.write("VASP=/home-yg/Soft/Vasp5.4/vasp_std\n")
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

                fout.write("cat >POSCAR<<EOF\n")
                self.poscar.to_poscar(fout)
                fout.write("EOF\n")
                fout.write("# scf\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_scf)
                fout.write("EOF\n")
                fout.write("cat > KPOINTS<<EOF\n")
                #self.kpoints.to_kpoints(fout)
                fout.write(kpoints_scf)
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_VASP_STD\n")

                fout.write("# nscf\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_nscf)
                fout.write("EOF\n")
                fout.write("cat >KPOINTS<<EOF\n")
                #self.kpoints.to_kpoints(fout)
                fout.write(kpoints_nscf)
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_VASP_STD\n")


                fout.write("# band structure\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_bands)
                fout.write("EOF\n")
                fout.write("cat >KPOINTS<<EOF\n")
                #self.kpoints.to_kpoints(fout)
                fout.write(kpoints_bands)
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_VASP_STD\n")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash static.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static", server=self.run_params["server"])