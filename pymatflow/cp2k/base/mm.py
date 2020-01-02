#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


from pymatflow.cp2k.base.mm_forcefield import cp2k_mm_forcefield
from pymatflow.cp2k.base.mm_neighbor_lists import cp2k_mm_neighbor_lists
from pymatflow.cp2k.base.mm_periodic_efield import cp2k_mm_periodic_efield
from pymatflow.cp2k.base.mm_poisson import cp2k_mm_poisson
from pymatflow.cp2k.base.mm_print import cp2k_mm_print

"""
usage:
"""

# ============================================
# CP2K / MM
#=============================================

class cp2k_mm:
    """

    """
    def __init__(self):
        """
        """
        self.params = {
                }
        self.status = False

        self.forcefield = cp2k_mm_force_field()
        self.neighbor_lists = cp2k_mm_neighbor_lists()
        self.periodic_efield = cp2k_mm_periodic_efield()
        self.poisson = cp2k_mm_poisson()
        self.printout = cp2k_mm_print()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&MM\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.forcefield.status == True:
            self.forcefield.to_input(fout)
        if self.neighbor_lists.status == True:
            self.neighbor_lists.to_input(fout)
        if self.periodic_efield.status == True:
            self.periodic_efield.to_input(fout)
        if self.poisson.status == True:
            self.poisson.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t&END DFT\n")
    
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "FORCEFIELD":
                self.forcefield.set_params({item: params[item]})
            elif item.split("-")[1] == "NEIGHBOR_LISTS":
                self.neighbor_lists.set_params({item: params[item]})
            elif item.split("-")[1] == "PERIODIC_EFIELD":
                self.periodic_efield.set_params({item: params[item]})
            elif item.split("-")[1] == "POISSON":
                self.poisson.set_params({item: params[item]})
            elif item.split("-")[1] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass
