#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import sys
import os
import shutil
import seekpath
import numpy as np
import pymatgen as mg

from pymatflow.cp2k.base.xyz import cp2k_xyz

from pymatflow.cp2k.base.glob import cp2k_glob
from pymatflow.cp2k.base.force_eval import cp2k_force_eval

"""
Usage:
    python phonon_cp2k.py xxx.xyz 
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.

Dependencies:
    pip install --user phonopy
    pip install --user cp2k_tools

Note:
    phonopy read the xxx.inp and it can only read the system structure
    by COORD specified in SUBSYS. So I can not use TOPOLOGY.
    PLUS: only scaled coordinates are currently supported!

    phonopy and seekpath both use spglib to decide the space group.

References:
    https://www.cp2k.org/exercises:2018_uzh_cmest:phonon_calculation
"""


class phonopy_run:
    """
    """
    def __init__(self, xyz_f):
        self.glob = cp2k_glob()
        self.force_eval = cp2k_force_eval(xyz_f)
        
        self.glob.basic_setting(run_type="ENERGY_FORCE")
        self.force_eval.basic_setting()

        self.supercell_n = [1, 1, 1]


    def phonopy(self, directory="tmp-cp2k-phonopy",
            mpi="", runopt="gen", force_eval={}):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            
            os.chdir(directory)
            shutil.copyfile("../%s" % self.force_eval.subsys.xyz.file, "%s" % self.force_eval.subsys.xyz.file)

            self.force_eval.set_params(force_eval)
            inp_name = "phonon.inp"
            with open(inp_name, 'w') as fout:
                self.glob.to_input(fout)
                #fout.write("\n")
                fout.write("&FORCE_EVAL\n")
                fout.write("\tMETHOD Quickstep\n")
            # subsys
            self.to_subsys_phonopy(inp_name)
            # end subsys
            with open(inp_name, 'a') as fout:
                # dft
                self.force_eval.dft.to_input(fout)
                # end dft
                fout.write("&END FORCE_EVAL\n")


            # build the phonopy running files 
            os.system("phonopy --cp2k -c %s -d --dim='%d %d %d'" % (inp_name, self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))
            # now phonon-supercell-00x.inp is generated which will be used to construct input for cp2k
            # in the past the above command will generate the supercell-00x.inp which can not be the
            # input file of cp2k directly, and need us to construct the final input for cp2k for every
            # displacement. but after an update of phonopy to version 2.4.2. the enerated phonon-supercell-00x.inp
            # can be run by cp2k directly! so we comment those old code for constructing the final input file.
            os.system("ls | grep 'phonon-supercell-' > geo.data")
            disps = []
            with open("geo.data", 'r') as fin:
                for line in fin:
                    disps.append(line.split(".")[0].split("-")[2])
            
            #for disp in disps:
            #    in_name = "phonon-supercell-%s.inp" % disp
            #    if os.path.exists(in_name) is not True:
            #        break
            #    tmp_file = "phonon-supercell-%s.tmp.txt" % disp
            #    shutil.copyfile(in_name, tmp_file)
            #    # important: different disp calculation should have different PROJECT name
            #    self.glob.params["PROJECT"] = "abinitio" + "-supercell-" + disp
            #    with open(in_name, 'w') as fout:
            #        self.glob.to_input(fout)
            #        fout.write("\n")
            #        fout.write("&FORCE_EVAL\n")
            #        fout.write("\tMETHOD Quickstep\n")
            #        fout.write("\t&SUBSYS\n")
            #    self.print_kinds(in_name)
            #    os.system("cat %s | sed '1d;2d;3d;4d;5d;6d;7d' | sed '$d' | sed '$d' | sed '$d' | sed '$d' | sed '$d' >> %s" % (tmp_file, in_name))
            #    with open(in_name, 'a') as fout:
            #        fout.write("\t&END SUBSYS\n")
            #        # dft
            #        self.force_eval.dft.to_input(fout)
            #        # end dft
            #        fout.write("\t&PRINT\n")
            #        fout.write("\t\t&FORCES\n")
            #        fout.write("\t\t\tFILENAME forces\n")
            #        fout.write("\t\t&END FORCES\n")
            #        fout.write("\t&END PRINT\n")
            #        fout.write("&END FORCE_EVAL\n")
            os.chdir("../")


            #
            # generate yhbatch file
            with open(os.path.join(directory, "phonopy-job.sub"), 'w') as fout:
                fout.write("#!/bin/bash\n\n")
                for disp in disps:
                    fout.write("yhrun -N 1 -n 24 cp2k.psmp -in phonon-supercell-%s.inp > phonon-supercell-%s.inp.out\n" % (disp, disp))

            # generate the result analysis bash script and necessary config files
            os.chdir(directory) 
            with open("mesh.conf", 'w') as fout:
                fout.write("ATOM_NAME =")
                for element in self.force_eval.subsys.xyz.specie_labels:
                    fout.write(" %s" % element)
                fout.write("\n")
                fout.write("DIM = %d %d %d\n" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))
                fout.write("MP = 8 8 8\n")
            with open("pdos.conf", 'w') as fout:
                fout.write("ATOM_NAME =")
                for element in self.force_eval.subsys.xyz.specie_labels:
                    fout.write(" %s" % element)
                fout.write("\n")
                fout.write("DIM = %d %d %d\n" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))
                fout.write("MP = 8 8 8\n")
                fout.write("PDOS = 1 2, 3 4 5 5\n")
            # plot the phonon band
            # 注意设置Primitive Axis要设置正确!
            with open("band.conf", 'w') as fout:
                fout.write("ATOM_NAME =")
                for element in self.force_eval.subsys.xyz.specie_labels:
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
                lattice = self.force_eval.subsys.xyz.cell #[self.force_eval.subsys.xyz.cell[0:3], self.force_eval.subsys.xyz.cell[3:6], self.force_eval.subsys.xyz.cell[6:9]]
                positions = []
                numbers = []
                a = np.sqrt(self.force_eval.subsys.xyz.cell[0][0]**2 + self.force_eval.subsys.xyz.cell[0][1]**2 + self.force_eval.subsys.xyz.cell[0][2]**2)
                b = np.sqrt(self.force_eval.subsys.xyz.cell[1][0]**2 + self.force_eval.subsys.xyz.cell[1][1]**2 + self.force_eval.subsys.xyz.cell[1][2]**2)
                c = np.sqrt(self.force_eval.subsys.xyz.cell[2][0]**2 + self.force_eval.subsys.xyz.cell[2][1]**2 + self.force_eval.subsys.xyz.cell[2][2]**2)
                for atom in self.force_eval.subsys.xyz.atoms:
                    positions.append([atom.x / a, atom.y / b, atom.z / c])
                    numbers.append(self.force_eval.subsys.xyz.specie_labels[atom.name])
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
                    fout.write(" $%s$" % point)
                point = kpoints_seekpath["path"][0][1]
                if point == "GAMMA":
                    fout.write(" $\Gamma$")
                else:
                    fout.write(" $%s$" % point)
                for i in range(1, len(kpoints_seekpath["path"])):
                    if kpoints_seekpath["path"][i][0] == kpoints_seekpath["path"][i-1][1]:
                        point = kpoints_seekpath["path"][i][1]
                        if point == "GAMMA":
                            fout.write(" $\Gamma$")
                        else:
                            fout.write(" $%s$" % point)
                    else:
                        point = kpoints_seekpath["path"][i][0]
                        if point == "GAMMA":
                            fout.write(" $\Gamma$")
                        else:
                            fout.write(" $%s$" % point)
                        point = kpoints_seekpath["path"][i][1]
                        if point == "GAMMA":
                            fout.write(" $\Gamma$")
                        else:
                            fout.write(" $%s$" % point)
                fout.write("\n")
            
            with open("phonopy-analysis.sh", 'w') as fout:
                fout.write("#!/bin/bash\n\n")
                fout.write("# get the FORCE_SETS\n")
                base_project_name = "ab-initio"
                phonopy_command = "phonopy --cp2k -f "
                for disp in disps:
                    # important: different disp calculation should have different PROJECT name
                    #f_name = "abinitio" + "-supercell-" + disp + "-forces-1_0.xyz"
                    f_name = "ab-initio" + "-supercell-" + disp + "-forces-1_0.xyz"
                    phonopy_command = phonopy_command + f_name + " "
                fout.write("%s\n" % phonopy_command)
                fout.write("# plot The density of states (DOS)\n")
                fout.write("phonopy --cp2k -p mesh.conf -c %s\n" % inp_name)
                fout.write("# Thermal properties are calculated with the sampling mesh by:\n")
                fout.write("phonopy --cp2k -t mesh.conf -c %s\n" % inp_name)
                fout.write("# Thermal properties can be plotted by:\n")
                fout.write("phonopy --cp2k -t -p mesh.conf -c %s\n" % inp_name)
                fout.write("# calculate Projected DOS and plot it\n")
                fout.write("phonopy --cp2k -p pdos.conf -c %s\n" % inp_name)
                fout.write("# get the band structure\n")
                fout.write("phonopy --cp2k -c %s -p band.conf\n" % inp_name)
            os.chdir("../")
            # end generate the result analysis bash script and necessray config file

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            disps = []
            with open("geo.data", 'r') as fin:
                for line in fin:
                    disps.append(line.split(".")[0].split("-")[2])
            for disp in disps:
                in_name = "supercell-%s.inp" % disp
                os.system("cp2k.psmp -in phonon-supercell-%s.inp | tee phonon-supercell-%s.inp.out" % (disp, disp))
            os.chdir("../")


    def to_subsys_phonopy(self, fname):
        cell = self.force_eval.subsys.xyz.cell
        with open(fname, 'a') as fout:
            fout.write("\t&SUBSYS\n")
            for element in self.force_eval.subsys.xyz.specie_labels:
                fout.write("\t\t&KIND %s\n" % element)
                fout.write("\t\t\tBASIS_SET DZVP-MOLOPT-SR-GTH\n")
                fout.write("\t\t\tPOTENTIAL GTH-PBE\n")
                fout.write("\t\t&END KIND\n")
            fout.write("\t\t&CELL\n")
            #fout.write("\t\t\tABC %f %f %f\n" % (cell[0], cell[4], cell[8]))
            # unit cell here can only be specified via ABC when doing phonopy calculation
            fout.write("\t\t\tABC %f %f %f\n" % (cell[0][0], cell[1][1], cell[2][2]))
            fout.write("\t\t&END CELL\n")
            #fout.write("\t\t&TOPOLOGY\n")
            #fout.write("\t\t\tCOORD_FILE_FORMAT xyz\n")
            #fout.write("\t\t\tCOORD_FILE_NAME %s\n" % sys.argv[1])
            #fout.write("\t\t&END TOPOLOGY\n")
            fout.write("\t\t&COORD\n")
            fout.write("\t\t\tSCALED .TRUE.\n")
            for atom in self.force_eval.subsys.xyz.atoms:
                #fout.write("\t\t\t%s\t%f\t%f\t%f\n" % (atom.name, atom.x/cell[0], atom.y/cell[4], atom.z/cell[8]))
                fout.write("\t\t\t%s\t%f\t%f\t%f\n" % (atom.name, atom.x/cell[0][0], atom.y/cell[1][1], atom.z/cell[2][2]))
            fout.write("\t\t&END COORD\n")
            fout.write("\t&END SUBSYS\n")
            fout.write("\n")

    def print_kinds(self, fname):
        with open(fname, 'a') as fout:
            for element in self.force_eval.subsys.xyz.specie_labels:
                fout.write("\t\t&KIND %s\n" % element)
                fout.write("\t\t\tBASIS_SET DZVP-MOLOPT-SR-GTH\n")
                fout.write("\t\t\tPOTENTIAL GTH-PBE\n")
                fout.write("\t\t&END KIND\n")

