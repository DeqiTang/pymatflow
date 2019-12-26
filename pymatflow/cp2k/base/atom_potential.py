#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from pymatflow.cp2k.base.atom_print import cp2k_atom_print

"""
Usage:
"""


class cp2k_atom_potential_ecp:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&ECP\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END ECP\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_potential_gth_potential:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&GTH_POTENTIAL\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END GTH_POTENTIAL\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_potential:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.ecp = cp2k_atom_potential_ecp()
        self.gth_potential = cp2k_atom_potential_gth_potential()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&POTENTIAL\n")
        for item in self.params:
            fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.ecp.status == True:
            self.ecp.to_input(fout)
        if self.gth_potential.status == True:
            self.gth_potential.to_input(fout)
        fout.write("\t&END POTENTIAL\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "ECP":
                self.ecp.set_params({item: params[item]})
            elif item.split("-")[1] == "GTH_POTENTIAL":
                self.gth_potential.set_params({item: params[item]})
            else:
                pass



