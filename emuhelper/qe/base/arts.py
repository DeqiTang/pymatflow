#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import pymatgen as mg
import os
import sys
import re

from emuhelper.base.xyz import base_xyz


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

        self.nks = 4

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
        fout.write("%f %f %f\n" % (cell[0], cell[1], cell[2]))
        fout.write("%f %f %f\n" % (cell[3], cell[4], cell[5]))
        fout.write("%f %f %f\n" % (cell[6], cell[7], cell[8]))
        fout.write("\n")
        fout.write("ATOMIC_POSITIONS angstrom\n")
        for atom in self.xyz.atoms:
            fout.write("%s\t%f\t%f\t%f\n" % (atom.name, atom.x, atom.y, atom.z))
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
        elif self.kpoints_option == "crystal_b":
            fout.write("K_POINTS %s\n" % self.kpoints_option)
            fout.write("%d\n" % self.nks)
            # Gamma-K-M-Gamma
            fout.write("0.0000000000     0.0000000000     0.0000000000     20\n")
            fout.write("0.3333333333     0.3333333333     0.0000000000     10\n")
            fout.write("0.0000000000     0.5000000000     0.0000000000     17\n")
            fout.write("0.0000000000     0.0000000000     0.0000000000     0\n")


    def set_kpoints(self, kpoints_mp=[1, 1, 1, 0, 0, 0], option="automatic"):
        """
        TODO: 
            considering using seekpath to get the kpoints automatically from structure
            https://github.com/giovannipizzi/seekpath/tree/develop/seekpath
        """
        if option == "automatic":
            self.kpoints_mp = kpoints_mp
        if option == "crystal_b":
            self.kpoints_option = option
            return
