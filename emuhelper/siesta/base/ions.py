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


class siesta_ions:
    """
    """
    def __init__(self):
        self.params = {
                "WriteCoorXmol": None,
                }
        self.md = {
                "TypeOfRun": None,
                "VariableCell": None,
                "ConstantVolume": None,
                "MaxForceTol": None,
                "MaxStressTol": None,
                "Steps": None,
                "MaxDispl": None,
                "PreconditionVariableCell": None,
                }


    def to_fdf(self, fout):
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item, str(self.params[item])))
        for item in self.md:
            if self.md[item] is not None:
                if item == "InitialTemperature":
                    fout.write("MD.%s %s K\n" % (item, str(self.md[item])))
                elif item == "TargetTemperature":
                    fout.write("MD.%s %s K\n" % (item, str(self.md[item])))
                elif item == "LengthTimeStep":
                    fout.write("MD.%s %s fs\n" % (item, str(self.md[item])))
                elif item == "MaxForceTol":
                    fout.write("MD.%s %s eV/Ang\n" % (item, str(self.md[item])))
                elif item == "MaxStressTol":
                    fout.write("MD.%s %s GPa\n" % (item, str(self.md[item])))
                elif item == "MaxDispl":
                    fout.write("MD.%s %s Bohr\n" % (item, str(self.md[item])))
                elif item == "PreconditionVariableCell":
                    fout.write("MD.%s %s Ang\n" % (item, str(self.md[item])))
                else:
                    fout.write("MD.%s %s\n" % (item, str(self.md[item])))
        #
        fout.write("\n")

    def basic_setting(self):
        self.md["TypeOfRun"] = "CG"   # CG, Broyden, 
        self.md["VariableCell"] = "false"
        self.md["ConstantVolume"] = "true"
        self.md["MaxForceTol"] = 0.001 # eV/Ang
        self.md["MaxStressTol"] = 0.01 # GPa
        self.md["Steps"] = 60
        self.md["MaxDispl"] = 0.2 # Bohr
        self.md["PreconditionVariableCell"] = 5 # Ang

        self.params["WriteCoorXmol"] = "true"
        self.params["WriteMDXmol"] = "true"
