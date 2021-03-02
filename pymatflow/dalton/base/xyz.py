#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

import pymatflow.base as base
from emuhelper.base.atom import Atom
from emuhelper.base.xyz import BaseXyz

"""
"""


class dalton_xyz(BaseXyz):
    """
    a representation of xyz file
    """
    def __init__(self, xyz_f):
        super().__init__(xyz_f)

    def to_dalton(self, fname):
        with open(fname, 'a') as fout:
            fout.write("Atomtypes=%d Angstrom\n" % self.nspecies)
            for element in self.specie_labels:
                charge = float(base.element[element].number)
                num_atoms = 0
                for atom in self.atoms:
                    if atom.name == element:
                        num_atoms += 1
                fout.write("Charge=%f Atoms=%d\n" % (charge, num_atoms))
                for atom in self.atoms:
                    if atom.name == element:
                        fout.write("%s %f %f %f\n" % (atom.name, atom.x, atom.y, atom.z))
            # end 

