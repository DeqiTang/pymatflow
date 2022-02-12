"""
Phonopy calc
"""
import numpy as np
import sys
import os
import shutil

from pymatflow.remote.server import server_handle

import pymatflow.base as base
from pymatflow.siesta.siesta import Siesta
#from pymatflow.siesta.base.system import SiestaSystem
#from pymatflow.siesta.base.electrons import SiestaElectrons



"""
Note:
    参考:
    https://atztogo.github.io/phonopy/siesta.html
"""

class PhonopyRun(Siesta):
    """
    """
    def __init__(self):
        super().__init__()
        #self.system = SiestaSystem()
        #self.electrons = SiestaElectrons()

        self.electrons.basic_setting()

        self.supercell_n = [1, 1, 1]


    def phonopy(self, directory="tmp-siesta-phonopy", inpname="phono-with-phonopy.fdf", output="phono-with-phonopy.out",
            runopt="gen", auto=0):
        """
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            shutil.copyfile(self.system.xyz.file, os.path.join(directory, os.path.basename(self.system.xyz.file)))
            for element in self.system.xyz.specie_labels:
                shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))


            # ok now we can use xyz class to extract information
            # from the xyz file: sys.argv[1]

            #xyz = siesta_xyz_phonopy()


            head_fdf_name = "head.fdf" # without sele.electrons now
            with open(os.path.join(directory, head_fdf_name), 'w') as fout:
                # we will add self.electrons.to_fdf after the supercell structure was append
                fout.write("SystemName %s\n" % self.system.name)
                fout.write("SystemLabel %s\n" % self.system.label)
                fout.write("NumberOfSpecies %s\n" % self.system.xyz.nspecies)

                fout.write(self.system.to_string())

                fout.write("%block ChemicalSpeciesLabel\n")
                for element in self.system.xyz.specie_labels:
                    fout.write("\t%d\t%d\t%s\n" % (self.system.xyz.specie_labels[element], base.element[element].number, element))
                fout.write("%endblock ChemicalSpeciesLabel\n")
                fout.write("\n")

                fout.write("%block PAO.BasisSizes\n")
                for element in self.system.xyz.specie_labels:
                    fout.write("\t%s\tDZP\n" % element)
                fout.write("%endblock PAO.BasisSizes\n")
                fout.write("\n")
                fout.write("# =========================================================\n")

            pos_fdf_name = "pos.fdf"
            with open(os.path.join(directory, pos_fdf_name), 'w') as fout:
                fout.write(self.system.to_string())

            # set up the Phonopy calculation
            os.chdir(directory)
            os.system("phonopy --siesta -d --dim='%d %d %d' -c %s" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2], pos_fdf_name))
            os.system("ls | grep 'supercell-' > pos.data")
            disp_dirs = []
            with open("pos.data", 'r') as fin:
                for line in fin:
                    disp_dirs.append(line.split(".")[0].split("-")[1])
            for disp in disp_dirs:
                os.mkdir("disp-%s" % disp)
                os.system("cat %s supercell-%s.fdf > ./disp-%s/supercell-%s.fdf" % (head_fdf_name, disp, disp, disp))
                # add electrons related setting to get the complete input for siesta
                with open("./disp-%s/supercell-%s.fdf" % (disp, disp), 'a') as fout:
                    fout.write("# =========================================================\n")
                    fout.write("\n\n")
                    fout.write(self.electrons.to_string())
                #
                os.system("cp *.psf ./disp-%s/" % disp)
                os.system("rm supercell-%s.fdf" % disp)
            os.chdir("../")


            # gen yhbatch script
            with open(os.path.join(directory, "phonopy-job.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                for disp in disp_dirs:
                    fout.write("cd disp-%s\n" % disp)
                    fout.write("yhrun $PMF_SIESTA < supercell-%s.fdf > supercell-%s.out\n" % (disp, disp))
                    fout.write("cd ../\n")

            # gen pbs script
            with open(os.path.join(directory, "phonopy-job.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for disp in disp_dirs:
                    fout.write("cd disp-%s\n" % disp)
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE siesta < supercell-%s.fdf > supercell-%s.out\n" % (disp, disp))
                    fout.write("cd ../\n")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            disp_dirs = []
            with open("pos.data", 'r') as fin:
                for line in fin:
                    disp_dirs.append(line.split(".")[0].split("-")[1])
            # run every disp
            for disp in disp_dirs:
                os.chdir("disp-%s" % disp)
                os.system("%s siesta < supercell-%s.fdf | tee supercell-%s.out" % (self.run_params["mpi"], disp, disp))
                os.chdir("../")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="phonopy-job", server=self.run_params["server"])
    #
    #
