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
"""


class siesta_transiesta:
    """
    """
    def __init__(self):
        self.ts = {
                "Voltage": None,
                "Elecs.Bulk": None,
                "Elecs.DM.Update": None,
                "Elecs.GF.ReUse": None,
                }

    def to_fdf(self, fout):
        for item in self.ts:
            if self.ts[item] is not None:
                fout.write("TS.%s %s\n" % (item, str(self.ts[item])))
 
        fout.write("%block TS.ChemPots\n")
        fout.write("Left\n")
        fout.write("Right\n")
        fout.write("%endblock TS.ChemPots\n")
 
        fout.write("%block TS.ChemPot.Left\n")
        fout.write("mu V/2\n")
        fout.write("contour.eq\n")
        fout.write("begin\n")
        fout.write("c-Left\n")
        fout.write("t-Left\n")
        fout.write("end\n")
        fout.write("%endblock TS.ChemPot.Left\n")

        fout.write("%block TS.ChemPot.Right\n")
        fout.write("mu -V/2\n")
        fout.write("contour.eq\n")
        fout.write("begin\n")
        fout.write("c-Right\n")
        fout.write("t-Right\n")
        fout.write("end\n")
        fout.write("%endblock TS.ChemPot.Right\n")
 

        fout.write("%block TS.Elecs\n")
        fout.write("Left\n")
        fout.write("Right\n")
        fout.write("%endblock TS.Elecs\n")
 
        fout.write("%block TS.Elec.Left\n")
        fout.write("HS siesta.TSHS\n")
        fout.write("chem-pot Left\n")
        fout.write("semi-inf-dir -a3\n")
        fout.write("elec-pos begin 1\n")
        fout.write("used-atoms 36\n")
        fout.write("%endblock TS.Elec.Left\n")
 
        fout.write("%block TS.Elec.Right\n")
        fout.write("HS siesta.TSHS\n")
        fout.write("chem-pot Right\n")
        fout.write("semi-inf-dir +a3\n")
        fout.write("elec-pos end -1\n")
        fout.write("used-atoms 36\n")
        fout.write("%endblock TS.Elec.Right\n")
 
        fout.write("TS.Contours.Eq.Pole    2.50000 eV\n")
        fout.write("%block TS.Contour.c-Left\n")
        fout.write("part circle\n")
        fout.write("from  -40.00000 eV + V/2 to -10. kT + V/2\n")
        fout.write("points 25\n")
        fout.write("method g-legendre\n")
        fout.write("%endblock TS.Contour.c-Left\n")
        fout.write("%block TS.Contour.t-Left\n")
        fout.write("part tail\n")
        fout.write("from prev to inf\n")
        fout.write("points 10\n")
        fout.write("method g-fermi\n")
        fout.write("%endblock TS.Contour.t-Left\n")
        fout.write("%block TS.Contour.c-Right\n")
        fout.write("part circle\n")
        fout.write("from  -40.00000 eV - V/2 to -10. kT - V/2\n")
        fout.write("points 25\n")
        fout.write("method g-legendre\n")
        fout.write("%endblock TS.Contour.c-Right\n")
        fout.write("%block TS.Contour.t-Right\n")
        fout.write("part tail\n")
        fout.write("from prev to inf\n")
        fout.write("points 10\n")
        fout.write("method g-fermi\n")
        fout.write("%endblock TS.Contour.t-Right\n")
 
        fout.write("TS.Elecs.Eta    0.0001000000 eV\n")
        fout.write("%block TS.Contours.nEq\n")
        fout.write("neq\n")
        fout.write("%endblock TS.Contours.nEq\n")
        fout.write("%block TS.Contour.nEq.neq\n")
        fout.write("part line\n")
        fout.write("from -|V|/2 - 5 kT to |V|/2 + 5 kT\n")
        fout.write("delta 0.01 eV\n")
        fout.write("method mid-rule\n")
        fout.write("%endblock TS.Contour.nEq.neq\n")

        fout.write("TS.WriteHS  .true.\n")
