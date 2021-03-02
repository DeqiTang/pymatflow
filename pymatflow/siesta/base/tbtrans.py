#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

from pymatflow.base.atom import Atom
from pymatflow.base.xyz import BaseXyz


"""
Usage:
"""


class SiestaTbtrans:
    """
    """
    def __init__(self):
        pass

    def to_fdf(self, fout):
        fout.write("TBT.Elecs.Eta    0.0001000000 eV\n")
        fout.write("%block TBT.Contours\n")
        fout.write("neq\n")
        fout.write("%endblock TBT.Contours\n")
 
        fout.write("%block TBT.Contour.neq\n")
        fout.write("part line\n")
        fout.write("from   -4.00000 eV to    4.00000 eV\n")
        fout.write("delta    0.04000 eV\n")
        fout.write("method mid-rule\n")
        fout.write("%endblock TBT.Contour.neq\n")
