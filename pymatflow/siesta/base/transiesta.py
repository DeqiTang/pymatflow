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


class SiestaTransiesta(SiestaVariableGroup):
    """
    """
    def __init__(self):
        super().__init__()

        self.set_param("TS.Voltage", 0, "eV")
        self.set_param("TS.Elecs.Bulk", "true")
        self.set_param("TS.Elecs.DM.Update", "cross-terms")
        self.set_param("TS.Elecs.GF.ReUse", "true")

    def to_string(self):
        out = ""
        out += "# Transiesta related parameters\n"
        out += super().to_string()
 
        out += "%block TS.ChemPots\n"
        out += "Left\n"
        out += "Right\n"
        out += "%endblock TS.ChemPots\n"
        out += "\n"
 
        out += "%block TS.ChemPot.Left\n"
        out += "mu V/2\n"
        out += "contour.eq\n"
        out += "begin\n"
        out += "c-Left\n"
        out += "t-Left\n"
        out += "end\n"
        out += "%endblock TS.ChemPot.Left\n"
        out += "\n"

        out += "%block TS.ChemPot.Right\n"
        out += "mu -V/2\n"
        out += "contour.eq\n"
        out += "begin\n"
        out += "c-Right\n"
        out += "t-Right\n"
        out += "end\n"
        out += "%endblock TS.ChemPot.Right\n"
        out += "\n"

        out += "%block TS.Elecs\n"
        out += "Left\n"
        out += "Right\n"
        out += "%endblock TS.Elecs\n"
        out += "\n"
 
        out += "%block TS.Elec.Left\n"
        out += "HS ../../electrodes/electrode-0/siesta.TSHS\n"
        out += "chem-pot Left\n"
        out += "semi-inf-dir -a3\n"
        out += "elec-pos begin 1\n"
        out += "used-atoms 36\n"
        out += "%endblock TS.Elec.Left\n"
        out += "\n"
 
        out += "%block TS.Elec.Right\n"
        out += "HS ../../electrodes/electrode-1/siesta.TSHS\n"
        out += "chem-pot Right\n"
        out += "semi-inf-dir +a3\n"
        out += "elec-pos end -1\n"
        out += "used-atoms 36\n"
        out += "%endblock TS.Elec.Right\n"
        out += "\n"
 
        out += "TS.Contours.Eq.Pole    2.50000 eV\n"
        out += "%block TS.Contour.c-Left\n"
        out += "part circle\n"
        out += "from  -40.00000 eV + V/2 to -10. kT + V/2\n"
        out += "points 25\n"
        out += "method g-legendre\n"
        out += "%endblock TS.Contour.c-Left\n"
        out += "%block TS.Contour.t-Left\n"
        out += "part tail\n"
        out += "from prev to inf\n"
        out += "points 10\n"
        out += "method g-fermi\n"
        out += "%endblock TS.Contour.t-Left\n"
        out += "%block TS.Contour.c-Right\n"
        out += "part circle\n"
        out += "from  -40.00000 eV - V/2 to -10. kT - V/2\n"
        out += "points 25\n"
        out += "method g-legendre\n"
        out += "%endblock TS.Contour.c-Right\n"
        out += "%block TS.Contour.t-Right\n"
        out += "part tail\n"
        out += "from prev to inf\n"
        out += "points 10\n"
        out += "method g-fermi\n"
        out += "%endblock TS.Contour.t-Right\n"
 
        out += "TS.Elecs.Eta    0.0001000000 eV\n"
        out += "%block TS.Contours.nEq\n"
        out += "neq\n"
        out += "%endblock TS.Contours.nEq\n"
        out += "%block TS.Contour.nEq.neq\n"
        out += "part line\n"
        out += "from -|V|/2 - 5 kT to |V|/2 + 5 kT\n"
        out += "delta 0.01 eV\n"
        out += "method mid-rule\n"
        out += "%endblock TS.Contour.nEq.neq\n"
        out += "\n"

        return out