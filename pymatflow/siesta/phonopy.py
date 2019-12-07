#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import seekpath
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
            os.system("phonopy --siesta -d --dim='%d %d %d' -c %s" % (self.supercelln[0], self.supercelln[1], self.supercelln[2], pos_fdf_name))
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
            
            with open("mesh.conf", 'w') as fout:
                fout.write("ATOM_NAME =")
                for element in self.system.xyz.specie_labels:
                    fout.write(" %s" % element)
                fout.write("\n")
                fout.write("DIM = %d %d %d\n" % (self.supercelln[0], self.supercelln[1], self.supercelln[2]))
                fout.write("MP = 8 8 8\n")
        
            # plot The density of states (DOS) 
            os.system("phonopy --siesta -p mesh.conf -c %s" % pos_fdf_name)
            # Thermal properties are calculated with the sampling mesh by:
            os.system("phonopy --siesta -t mesh.conf -c %s" % pos_fdf_name)
            # Thermal properties can be plotted by:
            os.system("phonopy --siesta -t -p mesh.conf -c %s" % pos_fdf_name)
        
            with open("pdos.conf", 'w') as fout:
                fout.write("ATOM_NAME =")
                for element in self.system.xyz.specie_labels:
                    fout.write(" %s" % element)
                fout.write("\n")
                fout.write("DIM = %d %d %d\n" % (self.supercelln[0], self.supercelln[1], self.supercelln[2]))
                fout.write("MP = 8 8 8\n")
                fout.write("PDOS = 1 2, 3 4 5 5\n")

            # calculate Projected DOS and plot it
            os.system("phonopy --siesta -p pdos.conf -c %s" % pos_fdf_name)
            # plot the phonon band
            with open("band.conf", 'w') as fout:
                fout.write("ATOM_NAME =")
                for element in self.system.xyz.specie_labels:
                    fout.write(" %s" % element)
                fout.write("\n")
                # the use of PRIMITIVE_AXES will find the primitive cell of the structure
                # and use it to analyse the phonon band structure
                # however, the use of primitive cell will not affect the q path setting
                # so whether we use PRIMITIVE cell or not, we can set the same q path
                fout.write("PRIMITIVE_AXES = AUTO\n") # we can also specify a matrix, but AUTO is recommended now in phonopy
                fout.write("GAMMA_CENTER = .TRUE.\n")
                fout.write("BAND_POINTS = 101\n")
                fout.write("BAND_CONNECTION = .TRUE.\n")

                fout.write("DIM = %d %d %d\n" % (self.supercelln[0], self.supercelln[1], self.supercelln[2]))
                #fout.write("BAND = 0.5 0.5 0.5 0.0 0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.0\n")
                fout.write("BAND =")
                # --------------
                # using seekpath to set q path
                # --------------
                lattice = [self.system.xyz.cell[0:3], self.system.xyz.cell[3:6], self.system.xyz.cell[6:9]]
                positions = []
                numbers = []
                a = np.sqrt(self.system.xyz.cell[0]**2 + self.system.xyz.cell[1]**2 + self.system.xyz.cell[2]**2)
                b = np.sqrt(self.system.xyz.cell[3]**2 + self.system.xyz.cell[4]**2 + self.system.xyz.cell[5]**2)
                c = np.sqrt(self.system.xyz.cell[6]**2 + self.system.xyz.cell[7]**2 + self.system.xyz.cell[8]**2)
                for atom in self.system.xyz.atoms:
                    positions.append([atom.x / a, atom.y / b, atom.z / c])
                    numbers.append(self.system.xyz.specie_labels[atom.name])
                structure = (lattice, positions, numbers)
                kpoints_seekpath = seekpath.get_path(structure)
                point = kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][0]]
                fout.write(" %f %f %f" % (point[0], point[1], point[2])) #self.arts..kpoints_seekpath["path"][0][0]
                point = kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][1]]
                fout.write(" %f %f %f" % (point[0], point[1], point[2])) #self.arts.kpoints_seekpath["path"][0][1]
                for i in range(1, len(kpoints_seekpath["path"])):
                    if kpoints_seekpath["path"][i][0] == kpoints_seekpath["path"][i-1][1]:
                        point = kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]]
                        fout.write(" %f %f %f" % (point[0], point[1], point[2])) #self.arts.kpoints_seekpath["path"][i][1]))
                    else:
                        point = kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][0]]
                        fout.write(" %f %f %f" % (point[0], point[1], point[2])) #self.arts.kpoints_seekpath["path"][i][0]))
                        point = kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]]
                        fout.write(" %f %f %f" % (point[0], point[1], point[2]))  #self.kpoints_seekpath["path"][i][1]))
                fout.write("\n")
                fout.write("BAND_LABELS =")
                point = kpoints_seekpath["path"][0][0]
                if point == "GAMMA":
                    fout.write(" $\Gamma$")
                else:
                    fout.write(" %s" % point)
                point = kpoints_seekpath["path"][0][1]
                if point == "GAMMA":
                    fout.write(" $\Gamma$")
                else:
                    fout.write(" %s" % point)
                for i in range(1, len(kpoints_seekpath["path"])):
                    if kpoints_seekpath["path"][i][0] == kpoints_seekpath["path"][i-1][1]:
                        point = kpoints_seekpath["path"][i][1]
                        if point == "GAMMA":
                            fout.write(" $\Gamma$")
                        else:
                            fout.write(" %s" % point)
                    else:
                        point = kpoints_seekpath["path"][i][0]
                        if point == "GAMMA":
                            fout.write(" $\Gamma$")
                        else:
                            fout.write(" %s" % point)
                        point = kpoints_seekpath["path"][i][1]
                        if point == "GAMMA":
                            fout.write(" $\Gamma$")
                        else:
                            fout.write(" %s" % point)
                fout.write("\n")
            os.system("phonopy --siesta -c %s -p band.conf" % pos_fdf_name)
            os.chdir("../")


    def gen_yh(self, inpname, output, directory="tmp-siesta-static", cmd="siesta"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % (cmd, inpname, output))

