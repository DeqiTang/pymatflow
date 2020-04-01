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
        self.magnetic_status = "spin-unpolarized" # "spin-polarized" "non-collinear"

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

            # gen llhpc script
            self.gen_llhpc(directory=directory, scriptname="static-scf.slurm", cmd="$PMF_VASP_STD")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="$PMF_VASP_STD", scriptname="static-scf.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # gen local bash script
            self.gen_bash(directory=directory, mpi=self.run_params["mpi"], cmd="$PMF_VASP_STD", scriptname="static-scf.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="$PMF_VASP_STD", scriptname="static-scf.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"])

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

            # gen llhpc script
            self.gen_llhpc(directory=directory, scriptname="static-nscf.slurm", cmd="$PMF_VASP_STD")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="$PMF_VASP_STD", scriptname="static-nscf.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # gen local bash script
            self.gen_bash(directory=directory, mpi=self.run_params["mpi"], cmd="$PMF_VASP_STD", scriptname="static-nscf.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="$PMF_VASP_STD", scriptname="static-nscf.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"])


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

            # gen llhpc script
            self.gen_llhpc(directory=directory, scriptname="static-bands.slurm", cmd="$PMF_VASP_STD")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="$PMF_VASP_STD", scriptname="static-bands.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # gen local bash script
            self.gen_bash(directory=directory, cmd="%s %PMF_VASP_STD" % self.run_params["mpi"], scriptname="static-bands.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="$PMF_VASP_STD", scriptname="static-bands.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"])

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
                    fout.write("yhrun $PMF_VASP_STD\n")
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
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi %s\n" % ("$PMF_VASP_STD"))
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
                    fout.write("yhrun $PMF_VASP_STD\n")
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
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi %s\n" % ("$PMF_VASP_STD"))
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
                    fout.write("yhrun %PMF_VASP_STD\n")
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
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi %s\n" % ("$PMF_VASP_STD"))
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

    def run(self, directory="tmp-vasp-static", runopt="gen", auto=0, kpath=None, kpath_intersections=15):
        """
        directory: a place for all the generated files

        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it
        Note: scf nscf(pdos, bands) in a single run
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.poscar.xyz.file, directory))

            self.incar.set_params({
                "IBRION": -1,
            })
            # scf
            #self.set_kpoints(option="automatic", kpoints_mp=kpoints_mp_scf)
            incar_scf = self.incar.to_string()
            kpoints_scf = self.kpoints.to_string()

            # nscf: pdos + bands
            self.incar.set_params({
                "ICHARG": 11,
                "LORBIT": 11,
                })            

            incar_nscf = self.incar.to_string()
            self.set_kpoints(option="bands", kpath=kpath, kpath_intersections=kpath_intersections)
            kpoints_nscf = self.kpoints.to_string()

            # gen llhpc script
            with open(os.path.join(directory, "static.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                fout.write("cat >POSCAR<<EOF\n")
                self.poscar.to_poscar(fout)
                fout.write("EOF\n")
                fout.write("# scf\n")
                fout.write("cat > INCAR<<EOF\n")
                fout.write(incar_scf)
                fout.write("EOF\n")
                fout.write("cat >KPOINTS<<EOF\n")
                #self.kpoints.to_kpoints(fout)
                fout.write(kpoints_scf)
                fout.write("EOF\n")

                if self.magnetic_status == "non-collinear":
                    fout.write("yhrun $PMF_VASP_NCL\n")
                else:
                    fout.write("yhrun $PMF_VASP_STD \n")
                fout.write("cp OUTCAR OUTCAR.scf\n")
                fout.write("cp vasprun.xml vasprun.xml.scf\n")

                fout.write("# nscf\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_nscf)
                fout.write("EOF\n")
                if self.incar.params["LHFCALC"] == ".TRUE." or self.incar.params["LHFCALC"] == "T" and float(self.incar.params["HFSCREEN"]) == 0.2:
                    fout.write("nk=`cat IBZKPT | head -n 2 | tail -n -1`\n")
                    nkpoint = 0
                    for i in range(len(kpath)-1):
                        if kpath[i][4] != "|":
                            nkpoint += kpath[i][4]
                    fout.write("nk=`echo \"${nk}+%d\" | bc`\n" % nkpoint)
                    fout.write("cat >KPOINTS<<EOF\n")
                    fout.write("Kpoint for HSE band structure\n")
                    fout.write("${nk}\n")
                    fout.write("EOF\n")
                    fout.write("cat IBZKPT | tail -n +3 >> KPOINTS\n")
                    fout.write("cat >> KPOINTS<<EOF\n")
                    #for kpoint in kpath:
                    #    fout.write("%f %f %f 0.0 !%s\n" % (kpoint[0], kpoint[1], kpoint[2], kpoint[3]))
                    for i in range(len(kpath)-1):
                        if kpath[i][4] != "|":
                            fout.write("%f %f %f 0.0 !%s\n" % (kpath[i][0], kpath[i][1], kpath[i][2], kpath[i][3]))
                            for j in range(kpath[i][4]-2):
                                x = (kpath[i+1][0] - kpath[i][0]) / (kpath[i][4]-1) * (j+1) + kpath[i][0]
                                y = (kpath[i+1][1] - kpath[i][1]) / (kpath[i][4]-1) * (j+1) + kpath[i][1]
                                z = (kpath[i+1][2] - kpath[i][2]) / (kpath[i][4]-1) * (j+1) + kpath[i][2]
                                fout.write("%f %f %f 0.0\n" % (x, y, z))
                            fout.write("%f %f %f 0.0 !%s\n" % (kpath[i+1][0], kpath[i+1][1], kpath[i+1][2], kpath[i+1][3]))
                        else:
                            continue
                    fout.write("EOF\n")
                else:
                    fout.write("cat >KPOINTS<<EOF\n")
                    fout.write(kpoints_nscf)
                    fout.write("EOF\n")
                
                if self.magnetic_status == "non-collinear":
                    fout.write("yhrun $PMF_VASP_NCL\n")
                else:
                    fout.write("yhrun $PMF_VASP_STD \n")
                fout.write("cp OUTCAR OUTCAR.nscf\n")
                fout.write('cp vasprun.xml vasprun.xml.nscf\n')


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
                if self.magnetic_status == "non-collinear":
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi $PMF_VASP_NCL \n")
                else:
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi $PMF_VASP_STD \n")
                fout.write("cp OUTCAR OUTCAR.scf\n")
                fout.write("cp vasprun.xml vasprun.xml.scf\n")

                fout.write("# nscf\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_nscf)
                fout.write("EOF\n")
                if self.incar.params["LHFCALC"] == ".TRUE." or self.incar.params["LHFCALC"] == "T" and float(self.incar.params["HFSCREEN"]) == 0.2:
                    fout.write("nk=`cat IBZKPT | head -n 2 | tail -n -1`\n")
                    nkpoint = 0
                    for i in range(len(kpath)-1):
                        if kpath[i][4] != "|":
                            nkpoint += kpath[i][4]
                    fout.write("nk=`echo \"${nk}+%d\" | bc`\n" % nkpoint)
                    fout.write("cat >KPOINTS<<EOF\n")
                    fout.write("Kpoint for HSE band structure\n")
                    fout.write("${nk}\n")
                    fout.write("EOF\n")
                    fout.write("cat IBZKPT | tail -n +3 >> KPOINTS\n")
                    fout.write("cat >> KPOINTS<<EOF\n")
                    #for kpoint in kpath:
                    #    fout.write("%f %f %f 0.0 !%s\n" % (kpoint[0], kpoint[1], kpoint[2], kpoint[3]))
                    for i in range(len(kpath)-1):
                        if kpath[i][4] != "|":
                            fout.write("%f %f %f 0.0 !%s\n" % (kpath[i][0], kpath[i][1], kpath[i][2], kpath[i][3]))                            
                            for j in range(kpath[i][4]-2):
                                x = (kpath[i+1][0] - kpath[i][0]) / (kpath[i][4]-1) * (j+1) + kpath[i][0]
                                y = (kpath[i+1][1] - kpath[i][1]) / (kpath[i][4]-1) * (j+1) + kpath[i][1]
                                z = (kpath[i+1][2] - kpath[i][2]) / (kpath[i][4]-1) * (j+1) + kpath[i][2]
                                fout.write("%f %f %f 0.0\n" % (x, y, z))
                            fout.write("%f %f %f 0.0 !%s\n" % (kpath[i+1][0], kpath[i+1][1], kpath[i+1][2], kpath[i+1][3]))
                        else:
                            continue
                    fout.write("EOF\n")          
                else:
                    fout.write("cat >KPOINTS<<EOF\n")
                    fout.write(kpoints_nscf)
                    fout.write("EOF\n")

                if self.magnetic_status == "non-collinear":
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi $PMF_VASP_NCL \n")
                else:
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi $PMF_VASP_STD \n")
                fout.write("cp OUTCAR OUTCAR.nscf\n")
                fout.write("cp vasprun.xml vasprun.xml.nscf\n")


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
                if self.magnetic_status == "non-collinear":
                    fout.write("%s $PMF_VASP_NCL \n" % self.run_params["mpi"])
                else:
                    fout.write("%s $PMF_VASP_STD \n" % self.run_params["mpi"])
                fout.write("cp OUTCAR OUTCAR.scf\n")
                fout.write("cp vasprun.xml vasprun.xml.scf\n")

                fout.write("# nscf\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_nscf)
                fout.write("EOF\n")
                if self.incar.params["LHFCALC"] == ".TRUE." or self.incar.params["LHFCALC"] == "T" and float(self.incar.params["HFSCREEN"]) == 0.2:
                    fout.write("nk=`cat IBZKPT | head -n 2 | tail -n -1`\n")
                    nkpoint = 0
                    for i in range(len(kpath)-1):
                        if kpath[i][4] != "|":
                            nkpoint += kpath[i][4]
                    fout.write("nk=`echo \"${nk}+%d\" | bc`\n" % nkpoint)
                    fout.write("cat >KPOINTS<<EOF\n")
                    fout.write("Kpoint for HSE band structure\n")
                    fout.write("${nk}\n")
                    fout.write("EOF\n")
                    fout.write("cat IBZKPT | tail -n +3 >> KPOINTS\n")
                    fout.write("cat >> KPOINTS<<EOF\n")
                    #for kpoint in kpath:
                    #    fout.write("%f %f %f 0.0 !%s\n" % (kpoint[0], kpoint[1], kpoint[2], kpoint[3]))
                    for i in range(len(kpath)-1):
                        if kpath[i][4] != "|":
                            fout.write("%f %f %f 0.0 !%s\n" % (kpath[i][0], kpath[i][1], kpath[i][2], kpath[i][3]))                            
                            for j in range(kpath[i][4]-2):
                                x = (kpath[i+1][0] - kpath[i][0]) / (kpath[i][4]-1) * (j+1) + kpath[i][0]
                                y = (kpath[i+1][1] - kpath[i][1]) / (kpath[i][4]-1) * (j+1) + kpath[i][1]
                                z = (kpath[i+1][2] - kpath[i][2]) / (kpath[i][4]-1) * (j+1) + kpath[i][2]
                                fout.write("%f %f %f 0.0\n" % (x, y, z))
                            fout.write("%f %f %f 0.0 !%s\n" % (kpath[i+1][0], kpath[i+1][1], kpath[i+1][2], kpath[i+1][3]))
                        else:
                            continue
                    fout.write("EOF\n")        
                else:
                    fout.write("cat >KPOINTS<<EOF\n")
                    fout.write(kpoints_nscf)
                    fout.write("EOF\n")

                if self.magnetic_status == "non-collinear":
                    fout.write("%s $PMF_VASP_NCL \n" % self.run_params["mpi"])
                else:
                    fout.write("%s $PMF_VASP_STD \n" % self.run_params["mpi"])
                fout.write("cp OUTCAR OUTCAR.nscf\n")
                fout.write("cp vasprun.xml vasprun.xml.nscf\n")


            # gen lsf_sz script
            with open(os.path.join(directory, "static.lsf_sz"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("APP_NAME=intelY_mid\n")
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
                if self.magnetic_status == "non-collinear":
                    fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_VASP_NCL\n")
                else:
                    fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_VASP_STD\n")
                fout.write("cp OUTCAR OUTCAR.scf\n")
                fout.write("cp vasprun.xml vasprun.xml.scf\n")
                fout.write("# nscf\n")
                fout.write("cat > INCAR<<EOF\n")
                #self.incar.to_incar(fout)
                fout.write(incar_nscf)
                fout.write("EOF\n")
                if self.incar.params["LHFCALC"] == ".TRUE." or self.incar.params["LHFCALC"] == "T" and float(self.incar.params["HFSCREEN"]) == 0.2:
                    fout.write("nk=`cat IBZKPT | head -n 2 | tail -n -1`\n")
                    nkpoint = 0
                    for i in range(len(kpath)-1):
                        if kpath[i][4] != "|":
                            nkpoint += kpath[i][4]
                    fout.write("nk=`echo \"${nk}+%d\" | bc`\n" % nkpoint)
                    fout.write("cat >KPOINTS<<EOF\n")
                    fout.write("Kpoint for HSE band structure\n")
                    fout.write("${nk}\n")
                    fout.write("EOF\n")
                    fout.write("cat IBZKPT | tail -n +3 >> KPOINTS\n")
                    fout.write("cat >> KPOINTS<<EOF\n")
                    #for kpoint in kpath:
                    #    fout.write("%f %f %f 0.0 !%s\n" % (kpoint[0], kpoint[1], kpoint[2], kpoint[3]))
                    for i in range(len(kpath)-1):
                        if kpath[i][4] != "|":
                            fout.write("%f %f %f 0.0 !%s\n" % (kpath[i][0], kpath[i][1], kpath[i][2], kpath[i][3]))                            
                            for j in range(kpath[i][4]-2):
                                x = (kpath[i+1][0] - kpath[i][0]) / (kpath[i][4]-1) * (j+1) + kpath[i][0]
                                y = (kpath[i+1][1] - kpath[i][1]) / (kpath[i][4]-1) * (j+1) + kpath[i][1]
                                z = (kpath[i+1][2] - kpath[i][2]) / (kpath[i][4]-1) * (j+1) + kpath[i][2]
                                fout.write("%f %f %f 0.0\n" % (x, y, z))
                            fout.write("%f %f %f 0.0 !%s\n" % (kpath[i+1][0], kpath[i+1][1], kpath[i+1][2], kpath[i+1][3]))
                        else:
                            continue
                    fout.write("EOF\n")           
                else:
                    fout.write("cat >KPOINTS<<EOF\n")
                    fout.write(kpoints_nscf)
                    fout.write("EOF\n")

                if self.magnetic_status == "non-collinear":
                    fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_VASP_NCL\n")
                else:
                    fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_VASP_STD\n")
                fout.write("cp OUTCAR OUTCAR.nscf\n")
                fout.write("cp vasprun.xml vasprun.xml.nscf\n")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash static.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static", server=self.run_params["server"])
