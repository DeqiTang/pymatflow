#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from pymatflow.base.atom import Atom
from pymatflow.base.xyz import base_xyz

"""
Usage:
    python converge_cutoff_cp2k.py xxx.xyz cutoff_min cutoff_max cutoff_step rel_cutoff
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.
"""


class cp2k_xyz(base_xyz):
    """
    a representation of xyz file
    """
    def __init__(self, xyz_f):
        super().__init__(xyz_f)


    def to_subsys(self, fname):
        cell = self.cell
        with open(fname, 'a') as fout:
            fout.write("\t&SUBSYS\n")
            for element in self.specie_labels:
                fout.write("\t\t&KIND %s\n" % element)
                fout.write("\t\t\tBASIS_SET DZVP-MOLOPT-SR-GTH\n")
                fout.write("\t\t\tPOTENTIAL GTH-PBE\n")
                fout.write("\t\t&END KIND\n")
            fout.write("\t\t&CELL\n")
            #fout.write("\t\t\tABC %f %f %f\n" % (cell[0], cell[4], cell[8]))
            fout.write("\t\t\tA %f %f %f\n" % (cell[0], cell[1], cell[2]))
            fout.write("\t\t\tB %f %f %f\n" % (cell[3], cell[4], cell[5]))
            fout.write("\t\t\tC %f %f %f\n" % (cell[6], cell[7], cell[8]))
            fout.write("\t\t&END CELL\n")
            fout.write("\t\t&TOPOLOGY\n")
            fout.write("\t\t\tCOORD_FILE_FORMAT xyz\n")
            fout.write("\t\t\tCOORD_FILE_NAME %s\n" % sys.argv[1])
            fout.write("\t\t&END TOPOLOGY\n")
            fout.write("\t&END SUBSYS\n")
            fout.write("\n")


