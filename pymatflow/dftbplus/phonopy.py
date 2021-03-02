#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

from pymatflow.remote.server import server_handle
from pymatflow.vasp.vasp import Vasp

"""
Usage:
    python phonon_with_phonopy_vasp.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and the POTCAR is in the directory.
    make sure the element in the POTCAR is in order of increasing the atom number of
    the element: 确保POTCAR中元素排列的顺序是原子序数从小到大
Note:
    参考: https://atztogo.github.io/phonopy/vasp.html

"""

class PhonopyRun(Vasp):
    """
    """
    def __init__(self):
        super().__init__()

        self.incar.set_runtype(runtype="static")

        self.supercell_n = [1, 1, 1]


    def phonopy(self, directory="tmp-vasp-phonopy", runopt="gen", auto=0):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.poscar.xyz.file, directory))

            # =======================================
            # Constructing the input file for VASP
            # =======================================


            # ======================================
            # Phonopy set up
            # ======================================

            # constructing the INCAR for the phonon calculation
            self.incar.params["ISMEAR"] = 0
            self.incar.params["SIGMA"] = 0.2
            self.incar.params["IBRION"] = -1
            self.incar.params["LWAVE"] = "F"
            self.incar.params["LCHARG"] = "F"

            #with open(os.path.join(directory, "INCAR"), 'w') as fout:
            #    self.incar.to_incar(fout)
            with open(os.path.join(directory, "KPOINTS"), "w") as fout:
                self.kpoints.to_kpoints(fout)

            ## Construct and run every POSCAR scf
            with open(os.path.join(directory, "POSCAR"), 'w') as fout:
                self.poscar.to_poscar(fout)

            os.chdir(directory)
            os.system("phonopy -d --dim='%d %d %d'" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))

            disps = self.get_disps("./")
            for disp in disps:
                os.mkdir("disp-%s" % (disp))
                os.chdir("disp-%s" % (disp))
                shutil.copyfile("../POSCAR-%s" % disp, "POSCAR")
                #shutil.copyfile("../INCAR", "INCAR")
                shutil.copyfile("../POTCAR", "POTCAR")
                shutil.copyfile("../KPOINTS", "KPOINTS")
                os.chdir("../")
            os.chdir("../") # end of input generation chdir outside of the directory

            # generate the llhpc script
            with open(os.path.join(directory, "phonopy-job.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                fout.write("cat > INCAR<<EOF\n")
                self.incar.to_incar(fout)
                fout.write("EOF\n")                
                for disp in disps:
                    fout.write("cd disp-%s\n" % disp)
                    fout.write("cp ../INCAR .\n")
                    fout.write("yhrun $PMF_VASP_STD\n")
                    fout.write("cd ../\n")

            # generate the pbs script
            with open(os.path.join(directory, "phonopy-job.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("cat > INCAR<<EOF\n")
                self.incar.to_incar(fout)
                fout.write("EOF\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for disp in disps:
                    fout.write("cd disp-%s\n" % disp)
                    fout.write("cp ../INCAR .\n")
                    #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi %s\n" % ("$PMF_VASP_STD"))
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s\n" % ("$PMF_VASP_STD"))
                    fout.write("cd ../\n")

            # generate the local bash script
            with open(os.path.join(directory, "phonopy-job.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n\n")
                fout.write("\n")
                fout.write("cat > INCAR<<EOF\n")
                self.incar.to_incar(fout)
                fout.write("EOF\n")
                for disp in disps:
                    fout.write("cd disp-%s\n" % disp)
                    fout.write("cp ../INCAR .\n")
                    fout.write("%s %s\n" % (self.run_params["mpi"], "$PMF_VASP_STD"))
                    fout.write("cd ../\n")

            # generate lsf_sz bash script
            with open(os.path.join(directory, "phonopy-job.lsf_sz"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("APP_NAME=%s\n" % self.run_params["queue"])
                fout.write("NP=%d\n" % self.run_params["nodes"]*self.run_params["ppn"]) #np)
                fout.write("NP_PER_NODE=%d\n" % self.run_params["ppn"]) #np_per_node)
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

                fout.write("cat > INCAR<<EOF\n")
                self.incar.to_incar(fout)
                fout.write("EOF\n")
                for disp in disps:
                    fout.write("cd disp-%s\n" % disp)
                    fout.write("cp ../INCAR .\n")
                    fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist $PMF_VASP_STD\n")
                    fout.write("cd ../\n")


            # generate lsf_sustc bash script
            with open(os.path.join(directory, "phonopy-job.lsf_sustc"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#BSUB -J %s\n" % self.run_params["jobname"])
                fout.write("#BSUB -q %s\n" % self.run_params["queue"])
                fout.write("#BSUB -n %s\n" % (self.run_params["nodes"] * self.run_params["ppn"])) #number of total cores
                fout.write("#BSUB -R \"span[ptile=%d]\"\n" % self.run_params["ppn"])
                fout.write("hostfile=`echo $LSB_DJOB_HOSTFILE`\n")
                fout.write("NP=`cat $hostfile | wc -l`\n")
                fout.write("cd $LS_SUBCWD\n")

                fout.write("cat > INCAR<<EOF\n")
                self.incar.to_incar(fout)
                fout.write("EOF\n")
                for disp in disps:
                    fout.write("cd disp-%s\n" % disp)
                    fout.write("cp ../INCAR .\n")
                    fout.write("mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP $PMF_VASP_STD\n")
                    fout.write("cd ../\n")

            # non-analytical term correction (optional)
            # 参见: https://atztogo.github.io/phonopy/vasp.html

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash phonopy-job.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="phonopy-job", server=self.run_params["server"])

    def get_disps(self, directory="./"):
        os.chdir(directory)
        os.system("ls | grep 'POSCAR-' > pos.data")

        disps = []
        with open("pos.data", 'r') as fin:
            for line in fin:
                disps.append(line.split("\n")[0].split("-")[1])
        return disps
