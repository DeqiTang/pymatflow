#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

from emuhelper.base.atom import Atom
from emuhelper.base.xyz import BaseXyz

"""
Usage:
"""

class orca_xyz(BaseXyz):
    """
    a representation of xyz file
    """
    def __init__(self, xyz_f):
        super().__init__(xyz_f)

    def to_orca(self, fname):
        with open(fname, 'a') as fout:
            for atom in self.atoms:
                fout.write("%s %f %f %f\n" % (atom.name, atom.x, atom.y, atom.z))

