#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

import pymatflow.base as base
from pymatflow.base.atom import Atom
from pymatflow.base.xyz import base_xyz

"""
Usage:

    make sure the xyz structure file and the pseudopotential file
    for all the elements of the system is in the directory.
"""


class abinit_xyz(base_xyz):
    """
    a representation of xyz file
    """
    def __init__(self, xyz_f):
        super().__init__(xyz_f)
    
    def to_abinit(self, fname):
        cell = self.cell
        with open(fname, 'a') as fout:
            fout.write("acell 1 1 1\n") # scaling with 1 means no actually scaling of rprim by acell
            fout.write("rprim\n")
            fout.write("%f %f %f\n" % (cell[0], cell[1], cell[2]))
            fout.write("%f %f %f\n" % (cell[3], cell[4], cell[5]))
            fout.write("%f %f %f\n" % (cell[6], cell[7], cell[8]))

            fout.write("ntypat %d\n" % self.nspecies)
            fout.write("natom %d\n" % self.natom)
            fout.write("typat\n")
            # abinit 不允许输入文件列数超过264, 因此如果原子数太多
            # 这里的typat要分多行列出
            # 利用余数设置如果一行超过30个原子就换行
            i = 0 
            for atom in self.atoms:
                fout.write("%d " % self.specie_labels[atom.name])
                if i % 30 == 29:
                    fout.write("\n")
                i += 1
            fout.write("\n")
            fout.write("znucl ")
            for element in self.specie_labels:
                fout.write(str(base.element[element].number))
                fout.write(" ")
            fout.write("\n")
            fout.write("\n")
            fout.write("xangst\n")
            for atom in self.atoms:
                fout.write("%f %f %f\n" % (atom.x, atom.y, atom.z))
            fout.write("\n")


