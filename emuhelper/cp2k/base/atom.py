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

class cp2k_atom_method_xc_wf_correlation:
    def __init__(self):
        self.params = {
                "METHOD": None,
                }
        self.params["METHOD"] = "RI_MP2_GPW"

    def to_xc(self, fout):
        fout.write("\t\t\t&WF_CORRELATION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&WFC_GPW\n")
        fout.write("\t\t\t\t\tCUTOFF 300\n")
        fout.write("\t\t\t\t\tREL_CUTOFF 50\n")
        fout.write("\t\t\t\t\tEPS_FILTER 1.0E-12\n")
        fout.write("\t\t\t\t\tEPS_GRID 1.0E-8\n")
        fout.write("\t\t\t\t&END WFC_GPW\n")
        fout.write("\t\t\t&END WF_CORRELATION\n")


class cp2k_atom_method_xc:
    def __init__(self):
        self.params = {
                }
        self.wf_correlation = cp2k_atom_method_xc_wf_correlation()

    def to_method(self, fout):
        fout.write("\t\t&XC\n")
        self.wf_correlation.to_xc(fout)
        fout.write("\t\t&END XC\n")


class cp2k_atom_method:
    def __init__(self):
        self.params = {
                }
        self.xc = cp2k_atom_method_xc()

    def to_atom(self, fout):
        fout.write("\t&METHOD\n")
        self.xc.to_method(fout)
        fout.write("\t&END METHOD\n")

class cp2k_atom:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.method = cp2k_atom_method()

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("&ATOM\n")
        for item in self.params:
            fout.write("\t%s %s\n" % (item, self.params[item]))
        self.method.to_atom(fout)
        fout.write("&END ATOM\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 1:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[0] == "METHOD":
                pass
