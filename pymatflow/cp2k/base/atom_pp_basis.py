#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

from pymatflow.cp2k.base.atom_print import cp2k_atom_print

"""
Usage:
"""



class cp2k_atom_pp_basis_basis:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&BASIS\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END BASIS\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_pp_basis:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.basis = cp2k_atom_pp_basis_basis()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&PP_BASIS\n")
        for item in self.params:
            fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.basis.status == True:
            self.basis.to_input(fout)
        fout.write("\t&END PP_BASIS\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "BASIS":
                self.basis.set_params({item: params[item]})
            else:
                pass


