#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

from pymatflow.base.atom import Atom
from pymatflow.base.xyz import BaseXyz

from pymatflow.siesta.group import SiestaVariableGroup

"""
Usage:
"""


class SiestaTbtrans(SiestaVariableGroup):
    """
    """
    def __init__(self):
        super().__init__()

    def to_string(self):
        out = ""
        out += "TBT.Elecs.Eta    0.0001000000 eV\n"
        out += "%block TBT.Contours\n"
        out += "neq\n"
        out += "%endblock TBT.Contours\n"
 
        out += "%block TBT.Contour.neq\n"
        out += "part line\n"
        out += "from   -4.00000 eV to    4.00000 eV\n"
        out += "delta    0.04000 eV\n"
        out += "method mid-rule\n"
        out += "%endblock TBT.Contour.neq\n"

        return out