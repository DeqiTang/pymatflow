#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil



"""
usage:
"""

# ============================================
# CP2K / PW_DFT / ITERATIVE_SOLVER
#=============================================

class cp2k_pw_dft_iterative_solver:
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
        fout.write("\t\t&ITERATIVE_SOLVER\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END ITERATIVE_SOLVER\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass
