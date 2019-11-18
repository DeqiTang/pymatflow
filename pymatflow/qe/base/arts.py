#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
import pymatgen as mg
import os
import sys
import re
import seekpath

from pymatflow.base.xyz import base_xyz


"""
usage:
"""

class qe_arts:
    """
    Control: &CELL, ATOMIC_SPECIES, ATOMIC_POSITIONS, 
             K_POINTS, CELL_PARAMETERS, CONSTRAINTS, 
             OCCUPATIONS, ATOMIC_FORCES
    """
    def __init__(self, xyz_f):
        self.xyz = base_xyz(xyz_f)

        self.cell_params = {
                "cell_dynamics": None,
                "press": None,
                "wmass": None,
                "cell_factor": None,
                "press_conv_thr": None,
                "cell_dofree": None,
                }

        self.kpoints_option = "automatic"
        self.kpoints_mp = [1, 1, 1, 0, 0, 0]
        self.kpoints_seekpath = None

        self.ifstatic = True # used to determine how to put atomic coordinates to input file

    def to_in(self, fout):
        # fout: a file stream for writing
        fout.write("&cell\n")
        fout.write("/\n")
        fout.write("\n")
        
        fout.write("ATOMIC_SPECIES\n")
        for element in self.xyz.specie_labels:
            tmp = os.listdir("./")
            pseudo_file = ""
            for f in tmp:
                match_string = "%s\." % element
                match = re.match(match_string, f)
                if match is not None and match.string.split(".")[-1] == 'UPF':
                    pseudo_file = match.string
                    break
            fout.write("%s %f %s\n" % (element, mg.Element(element).atomic_mass, pseudo_file))
        fout.write("\n")
        cell = self.xyz.cell
        fout.write("CELL_PARAMETERS angstrom\n")
        fout.write("%.9f %.9f %.9f\n" % (cell[0], cell[1], cell[2]))
        fout.write("%.9f %.9f %.9f\n" % (cell[3], cell[4], cell[5]))
        fout.write("%.9f %.9f %.9f\n" % (cell[6], cell[7], cell[8]))
        fout.write("\n")
        fout.write("ATOMIC_POSITIONS angstrom\n")
        if self.ifstatic == True:
            for atom in self.xyz.atoms:
                fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
        elif self.ifstatic == False:
            for atom in self.xyz.atoms:
                fout.write("%s\t%.9f\t%.9f\t%.9f" % (atom.name, atom.x, atom.y, atom.z))
                for fix in atom.fix:
                    if fix == True:
                        fout.write("\t0")
                    elif fix == False:
                        fout.write("\t1")
                fout.write("\n")
        else:
            print("===============================================\n")
            print("warning: qe.base.arts.to_in():\n")
            print("arts.ifstatic could only be True or False\n")
            sys.exit(1)
        fout.write("\n")
        
        # writing KPOINTS to the fout
        self.write_kpoints(fout)
        # =========================

    def write_kpoints(self, fout):
        # fout: a file stream for writing
        if self.kpoints_option == "automatic":
            fout.write("K_POINTS %s\n" % self.kpoints_option)
            fout.write("%d %d %d %d %d %d\n" % (
                self.kpoints_mp[0],
                self.kpoints_mp[1],
                self.kpoints_mp[2],
                self.kpoints_mp[3],
                self.kpoints_mp[4],
                self.kpoints_mp[5]
                ))
        elif self.kpoints_option == "gamma":
            fout.write("K_POINTS gamma\n")
        elif self.kpoints_option == "tpiba_b":
            fout.write("K_POINTS %s\n" % self.kpoints_option)
            nks = 2
            for i in range(1, len(self.kpoints_seekpath["path"])):
                if self.kpoints_seekpath["path"][i][0] == self.kpoints_seekpath["path"][i-1][1]:
                    nks = nks + 1
                else:
                    nks = nks + 2
            fout.write("%d\n" % nks)
            point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][0][0]]
            fout.write("%f %f %f %d  #%s\n" % (point[0], point[1], point[2], 5, self.kpoints_seekpath["path"][0][0]))
            point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][0][1]]
            fout.write("%f %f %f %d  #%s\n" % (point[0], point[1], point[2], 5, self.kpoints_seekpath["path"][0][1]))
            for i in range(1, len(self.kpoints_seekpath["path"])):
                if self.kpoints_seekpath["path"][i][0] == self.kpoints_seekpath["path"][i-1][1]:
                    point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][i][1]]
                    fout.write("%f %f %f %d  #%s\n" % (point[0], point[1], point[2], 5, self.kpoints_seekpath["path"][i][1]))
                else:
                    point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][i][0]]
                    fout.write("%f %f %f %d  #%s\n" % (point[0], point[1], point[2], 5, self.kpoints_seekpath["path"][i][0]))
                    point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][i][1]]
                    fout.write("%f %f %f %d  #%s\n" % (point[0], point[1], point[2], 5, self.kpoints_seekpath["path"][i][1]))
            #
        elif self.kpoints_option == "crystal_b":
            pass



    def set_kpoints(self, kpoints_mp=[1, 1, 1, 0, 0, 0], option="automatic"):
        """
        TODO: 
            considering using seekpath to get the kpoints automatically from structure
            https://github.com/giovannipizzi/seekpath/tree/develop/seekpath
        Note:
            "automatic" was controlled by kpoints_mp
            "gamma" was also handled internally
            other kpoints are set while seekpath

            seekpath will automatically find the primitive cell of the structure
            input by you. and the k points it gives is corresponding with that 
            primitive cell, if you stick to the use the original structure, you
            must know how to modify the k points to be applicable to your original
            structure.
        Plan:
            build a wrapper to the seekpath in a separate file[not decided now]
        """
        if option == "automatic":
            self.kpoints_option = option
            self.kpoints_mp = kpoints_mp
            return
        if option == "gamma":
            self.kpoints_option = option
            return
        # --------------
        # using seekpath
        # --------------
        lattice = [self.xyz.cell[0:3], self.xyz.cell[3:6], self.xyz.cell[6:9]]
        positions = []
        numbers = []
        a = np.sqrt(self.xyz.cell[0]**2 + self.xyz.cell[1]**2 + self.xyz.cell[2]**2)
        b = np.sqrt(self.xyz.cell[3]**2 + self.xyz.cell[4]**2 + self.xyz.cell[5]**2)
        c = np.sqrt(self.xyz.cell[6]**2 + self.xyz.cell[7]**2 + self.xyz.cell[8]**2)
        for atom in self.xyz.atoms:
            positions.append([atom.x / a, atom.y / b, atom.z / c])
            numbers.append(self.xyz.specie_labels[atom.name])
        structure = (lattice, positions, numbers)
        self.kpoints_seekpath = seekpath.get_path(structure)
        if option == "tpiba_b":
            self.kpoints_option = option
            return
    
    def basic_setting(self, ifstatic=True):
        self.ifstatic = ifstatic

    #

