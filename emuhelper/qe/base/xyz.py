#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import re
import pymatgen as mg

from emuhelper.base.atom import Atom
from emuhelper.base.xyz import base_xyz

"""
Usage:
    python geo_opt_not_vc_pwscf.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential
    file for all the elements in the system is in the directory.
"""


class qe_xyz(base_xyz):
    """
    a representation of xyz file
    """
    def __init__(self, xyz_f):
        super().__init__(xyz_f)


    def to_pwscf(self, fname):
        cell = self.cell
        with open(fname, 'a') as fout:
            fout.write("ATOMIC_SPECIES\n")
            for element in self.specie_labels:
                tmp = os.listdir("../")
                pseudo_file = ""
                for f in tmp:
                    match_string = "%s.*.UPF" % element
                    match = re.match(match_string, f)
                    if match is not None:
                        pseudo_file = match.string
                        break
                fout.write("%s %f %s\n" % (element, mg.Element(element).atomic_mass, pseudo_file))

            fout.write("CELL_PARAMETERS angstrom\n")
            fout.write("%f %f %f\n" % (cell[0], cell[1], cell[2]))
            fout.write("%f %f %f\n" % (cell[3], cell[4], cell[5]))
            fout.write("%f %f %f\n" % (cell[6], cell[7], cell[8]))
            fout.write("\n")
            fout.write("ATOMIC_POSITIONS angstrom\n")
            for atom in self.atoms:
                fout.write("%s\t%f\t%f\t%f\n" % (atom.name, atom.x, atom.y, atom.z))
            fout.write("\n")



