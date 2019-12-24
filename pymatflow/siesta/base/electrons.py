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
        self.kpoints_mp = [1, 1, 1]

    def to_fdf(self, fout):
        for item in self.params:
            if self.params[item] is not None:
                if item == "MeshCutoff":
                    fout.write("%s %s Ry\n" % (item, str(self.params[item])))
                elif item == "ElectronicTemperature":
                    fout.write("%s %s K\n" % (item, str(self.params[item])))
                else:
                    fout.write("%s %s\n" % (item, str(self.params[item])))
        for item in self.xc:
            if self.xc[item] is not None:
                fout.write("XC.%s %s\n" % (item, str(self.xc[item])))
        for item in self.dm:
            if self.dm[item] is not None:
                fout.write("DM.%s %s\n" % (item, str(self.dm[item])))
        fout.write("\n")
        fout.write("%block kgrid.MonkhorstPack\n")
        fout.write("%d 0 0 0.0\n" % self.kpoints_mp[0])
        fout.write("0 %d 0 0.0\n" % self.kpoints_mp[1])
        fout.write("0 0 %d 0.0\n" % self.kpoints_mp[2])
        fout.write("%endblock kgrid.MonkhorstPack\n")
        fout.write("\n")
        #
        fout.write("\n")

    def set_spin(self, spin="non-polarized"):
        self.params["Spin"] = spin

    def basic_setting(self):
        self.xc["functional"] = "GGA"
        self.xc["authors"] = "PBE"
        self.dm["Tolerance"] = 1.0-6
        self.dm["MixingWight"] = 0.1
        self.dm["NumberPulay"] = 8 # this can affect the convergence of scf
        self.dm["AllowExtrapolation"] = "true"
        self.dm["UseSaveDM"] = "false"
        self.params["SolutionMethod"] = "diagon"
        self.params["MeshCutoff"] = 300 #100
    
    def set_params(self, electrons):
    
        #for item in electrons:
        #    if item in self.params:
        #        self.params[item] = electrons[item]
        #        continue
        #    elif item.split(".")[1] in self.xc:
        #        self.xc[item.split(".")[1]] = electrons[item]
        #        continue
        #    elif item.split(".")[1] in self.dm:
        #        self.dm[item.split(".")[1]] = electrons[item]
        #        continue
        for item in electrons:
            if len(item.split(".")) == 1:
                self.params[item] = electrons[item]
            elif len(item.split(".")) == 2 and item.split(".")[0] == "XC":
                self.xc[item.split(".")[1]] = electrons[item]
            elif len(item.split(".")) == 2 and item.split(".")[0] == "DM":
                self.dm[item.split(".")[1]] = electrons[item]
        #

    
