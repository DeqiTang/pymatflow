#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import sys
import re
import os
import shutil
import seekpath
import numpy as np
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.qe.base.control import qe_control
from pymatflow.qe.base.system import qe_system
from pymatflow.qe.base.electrons import qe_electrons
from pymatflow.qe.base.arts import qe_arts



"""
Note:
    现在phonopy只支持设置ATOMIC_POSITIONS 为crystal类型
    参考:
    https://atztogo.github.io/phonopy/qe.html
"""




class phonopy_run:
    """
    Note:
        kpoints as well as energy cutoff both will have a significant
        influence on the precision of phonopy running.
        like when you have a small cell and use a small kpoints might
        result in a really useless phonon band. this is just as what
        we do in scf or opt running, generally a small structure requires
        large kpoint sets, while the larger one requires smaller kpoint
        sets.
    """
    def __init__(self, xyz_f):
        self.control = qe_control()
        self.system = qe_system()
        self.electrons = qe_electrons()
        #self.ions = qe_ions()
        self.arts = qe_arts(xyz_f)

        self.control.basic_setting("scf")
        self.system.basic_setting(self.arts)
        self.arts.basic_setting(ifstatic=True)

        self.supercell_n = [1, 1, 1]
        
    def phonopy(self, directory="tmp-qe-phonopy", pos_inpname="pos.in", head_inpname="head.in",
            mpi="", runopt="gen", control={}, system={}, electrons={}, 
            kpoints_option="automatic", kpoints_mp=[1, 1, 1, 0, 0, 0], supercell_n=[1, 1, 1]):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.UPF %s/" % directory)
            os.system("cp %s %s/" % (self.arts.xyz.file, directory))
           
            self.supercell_n = supercell_n


            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.set_occupations() which uses
            # self.system.set_occupations() to set them, as self.system.set_params() 
            # is suppressed from setting occupations related parameters
            self.set_occupations(system)
            self.control.set_params(control)
            self.system.set_params(system)
            self.electrons.set_params(electrons)
            self.arts.set_kpoints(option=kpoints_option, kpoints_mp=kpoints_mp)
            
            # must print print out forces and stress after scf
            # so that phonopy can parse the scf output file and
            # construct the FORCE CONSTANT MATIRX
            self.control.params["tprnfor"] = True
            self.control.params["tstress"] = True
            
                        
            with open(os.path.join(directory, head_inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.write_kpoints(fout)

            # set up the Phonopy calculation
            os.chdir(directory)
            os.system("cat %s > %s" % (head_inpname, pos_inpname))
            with open(pos_inpname, 'a') as fout:
                self.arts.to_in(fout, coordtype="crystal")
            os.system("phonopy --qe -d --dim='%d %d %d' -c %s" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2], pos_inpname))
            os.system("ls | grep 'supercell-' > pos.data")
            disp_dirs = []
            with open("pos.data", 'r') as fin:
                for line in fin:
                    disp_dirs.append(line.split(".")[0].split("-")[1])
            # IMPORTANT:
            # we need to overwrite the head_inpname file, reasons are bellow
            # here we must pay attention to system.params["nat"]
            # as system.params["nat"] is the value of the original structure
            # but after phonopy build the supercell the number of atoms changes
            # accordingly. so we must change the value of system.params["nat"]
            # according to self.supercell_n
            self.system.params["nat"] = self.system.params["nat"] * self.supercell_n[0] * self.supercell_n[1] * self.supercell_n[2]
            with open(head_inpname, 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.write_kpoints(fout)
            #
            for disp in disp_dirs:
                os.system("cat %s supercell-%s.in > supercell-%s-full.in" % (head_inpname, disp, disp))
                os.system("rm supercell-%s.in" % disp)
            os.chdir("../")
            # end build the phonopy

            # gen yhbatch script
            #self.gen_yh(directory=directory, inpname=inpname, output=output)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("%s pw.x < %s | tee %s" % (mpi, inpname, output))
            # run the dft
            for disp in disp_dirs:
                os.system("pw.x < supercell-%s-full.in | tee supercell-%s.out" % (disp, disp))

            # analyse the result
            import matplotlib.pyplot as plt

            os.system("phonopy --qe -f supercell-{001..%s}.out" % (disp_dirs[-1]))
            
            with open("mesh.conf", 'w') as fout:
                fout.write("ATOM_NAME =")
                for element in self.arts.xyz.specie_labels:
                    fout.write(" %s" % element)
                fout.write("\n")
                fout.write("DIM = %d %d %d\n" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))
                fout.write("MP = 8 8 8\n")
        
            # plot The density of states (DOS) 
            os.system("phonopy --qe -p mesh.conf -c %s" % pos_inpname)
            # Thermal properties are calculated with the sampling mesh by:
            os.system("phonopy --qe -t mesh.conf -c %s" % pos_inpname)
            # Thermal properties can be plotted by:
            os.system("phonopy --qe -t -p mesh.conf -c %s" % pos_inpname)
        
            with open("pdos.conf", 'w') as fout:
                fout.write("ATOM_NAME =")
                for element in self.arts.xyz.specie_labels:
                    fout.write(" %s" % element)
                fout.write("\n")
                fout.write("DIM = %d %d %d\n" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))
                fout.write("MP = 8 8 8\n")
                fout.write("PDOS = 1 2, 3 4 5 5\n")

            # calculate Projected DOS and plot it
            os.system("phonopy --siesta -p pdos.conf -c %s" % pos_inpname)
            # plot band structure
            with open("band.conf", 'w') as fout:
                fout.write("ATOM_NAME =")
                for element in self.arts.xyz.specie_labels:
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
                
                fout.write("DIM = %d %d %d\n" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))
                #fout.write("BAND = 0.5 0.5 0.5 0.0 0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.0\n")
                fout.write("BAND =")
                # --------------
                # using seekpath to set q path
                # --------------
                lattice = [self.arts.xyz.cell[0:3], self.arts.xyz.cell[3:6], self.arts.xyz.cell[6:9]]
                positions = []
                numbers = []
                a = np.sqrt(self.arts.xyz.cell[0]**2 + self.arts.xyz.cell[1]**2 + self.arts.xyz.cell[2]**2)
                b = np.sqrt(self.arts.xyz.cell[3]**2 + self.arts.xyz.cell[4]**2 + self.arts.xyz.cell[5]**2)
                c = np.sqrt(self.arts.xyz.cell[6]**2 + self.arts.xyz.cell[7]**2 + self.arts.xyz.cell[8]**2)
                for atom in self.arts.xyz.atoms:
                    positions.append([atom.x / a, atom.y / b, atom.z / c])
                    numbers.append(self.arts.xyz.specie_labels[atom.name])
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
                #
            #
            os.system("phonopy --qe -c %s -p band.conf" % pos_inpname)

            os.chdir("../")
    

    def set_occupations(self, system):
        """
            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.system.set_occupations() to 
            # set them, as self.system.set_params() is suppressed from setting
            # occupations related parameters
            # if occupations == None, use default smearing occupation. and 
            # if occupations == "tetrahedra" the value set for smearing and degauss is ignored.
            # if occupations == "smearing", the value of smearing and degauss
            # should be legacy, not None or other illegal values.
        """
        if "occupations" in system:
            if system["occupations"] == None: # user default setting of set_occupations()
                self.system.set_occupations()
            elif system["occupations"] == "tetrahedra":
                self.system.set_occupations(occupations="tetrahedra")
            elif system["occupations"] == "smearing":
                if "smearing" in system and "degauss" in system:
                    self.system.set_occupations(occupations="smearing", smearing=system["smearing"], degauss=system["degauss"])
                elif "smearing" in system:
                    self.system.set_occupations(occupations="smearing", smearing=system["smearing"])
                elif "degauss" in system:
                    self.system.set_occupations(occupations="smearing", degauss=system["degauss"])
                else:
                    self.system.set_occupations(occupations="smearing")
            elif system["occupations"] == "tetrahedra_lin":
                self.system.set_occupations(occupations="tetrahedra_lin")
            elif system["occupations"] == "tetrahedra_opt":
                self.system.set_occupations(occupations="tetrahedra_opt")
            elif system["occupations"] == "fixed":
                self.system.set_occupations(occupations="fixed")
            elif system["occupations"] == "from_input":
                self.system.set_occupations(occupations="from_input")

    def gen_yh(self, directory, inpname, output):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 pw.x < %s > %s\n" % (inpname, output))

