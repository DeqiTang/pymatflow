"""
in control of the electrons step related parameters
"""
import numpy as np
import sys
import os
import shutil

from pymatflow.siesta. group import SiestaVariableGroup

"""
Usage:
"""


class SiestaElectrons(SiestaVariableGroup):
    """
    """
    def __init__(self):
        super().__init__()
        self.incharge = [
                "MaxSCFIterations", "SolutionMethod", "MeshCutoff", "XC.functional", "XC.author",
                "DM.Tolerance", "DM.MixingWeight", "DM.NumberPulay", "DM.AllowExtrapolation",
                "DM.UseSaveDM",
                ]
        self.kpoints_mp = [1, 1, 1]

    def to_string(self):
        out = ""
        out += super().to_string()
        out += "\n"
        out += "\n"
        out += "%block kgrid.MonkhorstPack\n"
        out += "%d 0 0 0.0\n" % self.kpoints_mp[0]
        out += "0 %d 0 0.0\n" % self.kpoints_mp[1]
        out += "0 0 %d 0.0\n" % self.kpoints_mp[2]
        out += "%endblock kgrid.MonkhorstPack\n"
        out += "\n"
        #
        out += "\n"
        return out

    def set_spin(self, spin="non-polarized"):
        self.set_param("Spin", spin)

    def basic_setting(self):
        self.set_param("XC.functional", "GGA")
        self.set_param("XC.authors", "PBE")
        self.set_param("DM.Tolerance", 1.0-6)
        self.set_param("DM.MixingWight", 0.1)
        self.set_param("DM.NumberPulay", 8) # this can affect the convergence of scf
        self.set_param("DM.AllowExtrapolation", "true")
        self.set_param("DM.UseSaveDM", "false")
        self.set_param("SolutionMethod", "diagon")
        self.set_param("MeshCutoff", 300, "Ry") #100
        #

    
