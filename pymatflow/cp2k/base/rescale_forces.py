#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


"""
usage:
"""

# ============================================
# CP2K / RESCALE_FORCES
#=============================================

class cp2k_rescale_forces:
    """

    """
    def __init__(self):
        """
        """
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&RESCALE_FORCES\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END RESCALE_FORCES\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "ALMO_SCF":
                self.almo_scf.set_params({item: params[item]})
            else:
                pass
