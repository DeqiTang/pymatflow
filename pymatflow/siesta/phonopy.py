#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import seekpath
import pymatgen as mg


from pymatflow.siesta.siesta import siesta
#from pymatflow.siesta.base.system import siesta_system
#from pymatflow.siesta.base.electrons import siesta_electrons



"""
Note:
    参考:
    https://atztogo.github.io/phonopy/siesta.html
"""

class phonopy_run(siesta):
    """
    """
    def __init__(self):
        super().__init__()
        #self.system = siesta_system()
        #self.electrons = siesta_electrons()

        self.electrons.basic_setting()
            
        self.supercell_n = [1, 1, 1] 


    def phonopy(self, directory="tmp-siesta-phonopy", inpname="phono-with-phonopy.fdf", output="phono-with-phonopy.out",
            mpi="", runopt="gen",
            jobname="siesta-phonopy", nodes=1, ppn=32):
        """
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
        
            for element in self.system.xyz.specie_labels:
                shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))

            
            # ok now we can use xyz class to extract information 
            # from the xyz file: sys.argv[1]

            #xyz = siesta_xyz_phonopy()


            head_fdf_name = "head.fdf" # without sele.electrons now
            with open(os.path.join(directory, head_fdf_name), 'w') as fout:
                # we will add self.electrons.to_fdf after the supercell structure was append
                #self.electrons.to_fdf(fout)
                fout.write("SystemName %s\n" % self.system.name)
                fout.write("SystemLabel %s\n" % self.system.label)
                fout.write("NumberOfSpecies %s\n" % self.system.xyz.nspecies)

        
                import pymatgen as mg
                fout.write("%block ChemicalSpeciesLabel\n")
                for element in self.system.xyz.specie_labels:
                    fout.write("\t%d\t%d\t%s\n" % (self.system.xyz.specie_labels[element], mg.Element(element).number, element))
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
                self.system.to_fdf(fout)
            
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
                    self.electrons.to_fdf(fout)
                #
                os.system("cp *.psf ./disp-%s/" % disp)
                os.system("rm supercell-%s.fdf" % disp)
            os.chdir("../")


            # gen yhbatch script
            with open(os.path.join(directory, "phonopy-job.sub"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                for disp in disp_dirs:
                    fout.write("cd disp-%s\n" % disp)
                    fout.write("yhrun -N 1 -n 24 siesta < supercell-%s.fdf > supercell-%s.out\n" % (disp, disp))
                    fout.write("cd ../\n")

            # gen pbs script
            with open(os.path.join(directory, "phonopy-job.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % jobname)
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
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
                os.system("siesta < supercell-%s.fdf | tee supercell-%s.out" % (disp, disp))
                os.chdir("../")
            os.chdir("../")
    #
    #
