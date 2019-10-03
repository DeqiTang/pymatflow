#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg



"""
Usage:
"""


class siesta_electrons:
    """
    """
    def __init__(self):
        self.params = {
                "MaxSCFIterations": None,
                "SolutionMethod": None,
                "MeshCutoff": None,
                }
        self.xc = {
                "functional": None,
                "author": None,
                }
        self.dm = {
                "Tolerance": None,
                "MixingWeight": None,
                "NumberPulay": None,
                "AllowExtrapolation": None,
                "UseSaveDM": None,
                }


    def to_fdf(self, fout):
        for item in self.params:
            if self.params[item] is not None:
                if item == "MeshCutoff":
                    fout.write("%s %s Ry\n" % (item, str(self.params[item])))
                else:
                    fout.write("%s %s\n" % (item, str(self.params[item])))
        for item in self.xc:
            if self.xc[item] is not None:
                fout.write("XC.%s %s\n" % (item, str(self.xc[item])))
        for item in self.dm:
            if self.dm[item] is not None:
                fout.write("DM.%s %s\n" % (item, str(self.dm[item])))
        #
        fout.write("\n")
