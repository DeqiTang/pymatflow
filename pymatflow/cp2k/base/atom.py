#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


from pymatflow.cp2k.base.atom_ae_basis import cp2k_atom_ae_basis
from pymatflow.cp2k.base.atom_method import cp2k_atom_method
from pymatflow.cp2k.base.atom_optimization import cp2k_atom_optimization
from pymatflow.cp2k.base.atom_potential import cp2k_atom_potential
from pymatflow.cp2k.base.atom_powell import cp2k_atom_powell
from pymatflow.cp2k.base.atom_pp_basis import cp2k_atom_pp_basis
from pymatflow.cp2k.base.atom_print import cp2k_atom_print


"""
Usage:
"""


class cp2k_atom:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false
        
        self.ae_basis = cp2k_atom_ae_basis()
        self.method = cp2k_atom_method()
        self.optimization = cp2k_atom_optimization()
        self.potential = cp2k_atom_potential()
        self.powell = cp2k_atom_powell()
        self.pp_basis = cp2k_atom_pp_basis()
        self.printout = cp2k_atom_print()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("&ATOM\n")
        for item in self.params:
            fout.write("\t%s %s\n" % (item, self.params[item]))
        if self.ae_basis.status == True:
            self.ae_basis.to_input(fout)
        if self.method.status == True:
            self.method.to_input(fout)
        if self.optimization.status == True:
            self.optimization.to_input(fout)
        if self.potential.status == True:
            self.potential.to_input(fout)
        if self.powell.status == True:
            self.powell.to_input(fout)
        if self.pp_basis.status == True:
            self.pp_basis.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("&END ATOM\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 1:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[0] == "AE_BASIS":
                self.ae_basis.set_params({item: params[item]})
            elif item.split("-")[0] == "METHOD":
                self.method.set_params({item: params[item]})
            elif item.split("-")[0] == "OPTIMIZATION":
                self.optimization.set_params({item: params[item]})
            elif item.split("-")[0] == "POTENTIAL":
                self.potential.set_params({item: params[item]})
            elif item.split("-")[0] == "POWELL":
                self.powell.set_params({item: params[item]})
            elif item.split("-")[0] == "PP_BASIS":
                self.pp_basis.set_params({item: params[item]})
            elif item.split("-")[0] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass
