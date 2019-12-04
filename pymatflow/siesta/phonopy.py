#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


from pymatflow.siesta.base.system import siesta_system
from pymatflow.siesta.base.electrons import siesta_electrons



"""
Note:
    参考:
    https://atztogo.github.io/phonopy/siesta.html
"""

class phonopy_run:
    """
    """
    def __init__(self, xyz_f):
        self.system = siesta_system(xyz_f)
        self.electrons = siesta_electrons()

        self.electrons.basic_setting()
            
        self.supercelln = [1, 1, 1] 

    def phonopy(self, directory="tmp-siesta-phonopy", inpname="phono-with-phonopy.fdf", output="phono-with-phonopy.out",
            mpi="", runopt="gen", electrons={}, ions={}, kpoints_mp=[1, 1, 1], supercelln=[1, 1, 1]):
        """
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
        
            for element in self.system.xyz.specie_labels:
                shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))

            self.supercelln = supercelln
            self.electrons.kpoints_mp = kpoints_mp
            self.electrons.set_params(electrons)
            
            # OK now we can use XYZ class to extract information 
            # from the xyz file: sys.argv[1]

            #xyz = siesta_xyz_phonopy()


            head_fdf_name = "head.fdf"
            with open(os.path.join(directory, head_fdf_name), 'w') as fout:
                #self.system.to_fdf(fout)
                self.electrons.to_fdf(fout)
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

            pos_fdf_name = "pos.fdf"
            with open(os.path.join(directory, pos_fdf_name), 'w') as fout:
                self.system.to_fdf(fout)
            
            # set up the Phonopy calculation
            os.chdir(directory)
            os.system("phonopy --siesta -d --dim='%d %d %d' -c %s" % (self.supercelln[0], self.supercelln[1], self.supercelln[2], pos_fdf_name))
            os.system("ls | grep 'supercell-' > pos.data")
            disp_dirs = []
            with open("pos.data", 'r') as fin:
                for line in fin:
                    disp_dirs.append(line.split(".")[0].split("-")[1])
            for disp in disp_dirs:
                os.mkdir("disp-%s" % disp)
                os.system("cat %s supercell-%s.fdf > ./disp-%s/supercell-%s.fdf" % (head_fdf_name, disp, disp, disp))
                os.system("cp *.psf ./disp-%s/" % disp)
                os.system("rm supercell-%s.fdf" % disp)
            os.chdir("../")


            # gen yhbatch script
            #self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="siesta")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            # run every disp
            for disp in disp_dirs:
                os.chdir("disp-%s" % disp)
                os.system("siesta < supercell-%s.fdf | tee supercell-%s.out" % (disp, disp))
                os.chdir("../")

            # analyse the result
            import matplotlib.pyplot as plt
            # create FORCE_SETS
            os.system("phonopy --siesta -f disp-{001..%s}/%s.FA" % (disp_dirs[-1], "siesta"))
            # plot the phonon band
            with open("band.conf", 'w') as fout:
                fout.write("ATOM_NAME =")
                for element in self.system.xyz.specie_labels:
                    fout.write(" %s" % element)
                fout.write("\n")
                fout.write("DIM = %d %d %d\n" % (self.supercelln[0], self.supercelln[1], self.supercelln[2]))
                fout.write("BAND = 0.5 0.5 0.5 0.0 0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.0\n")
            os.system("phonopy --siesta -c %s -p band.conf" % pos_fdf_name)
            os.chdir("../")


    def gen_yh(self, inpname, output, directory="tmp-siesta-static", cmd="siesta"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % (cmd, inpname, output))

