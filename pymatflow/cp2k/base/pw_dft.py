#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil


from pymatflow.cp2k.base.pw_dft_control import cp2k_pw_dft_control
from pymatflow.cp2k.base.pw_dft_iterative_solver import cp2k_pw_dft_iterative_solver
from pymatflow.cp2k.base.pw_dft_mixer import cp2k_pw_dft_mixer
from pymatflow.cp2k.base.pw_dft_parameters import cp2k_pw_dft_parameters

"""
usage:
"""

# ============================================
# CP2K / PW_DFT
#=============================================

class cp2k_pw_dft:
    """
    """
    def __init__(self):
        """
        """
        self.params = {
                }
        self.status = False

        self.control = cp2k_pw_dft_control()
        self.iterative_solver = cp2k_pw_dft_iterative_solver()
        self.mixer = cp2k_pw_dft_mixer()
        self.parameters = cp2k_pw_dft_parameters()

        # basic setting
        self.control.status = True
        self.iterative_solver.status = True
        self.mixer.status = True
        self.parameters.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&PW_DFT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.control.status == True:
            self.control.to_input(fout)
        if self.iterative_solver.status == True:
            self.iterative_solver.to_input(fout)
        if self.mixer.status == True:
            self.mixer.to_input(fout)
        if self.parameters.status == True:
            self.parameters.to_input(fout)
        fout.write("\t&END PW_DFT\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "CONTROL":
                self.control.set_params({item: params[item]})
            elif item.split("-")[1] == "ITERATIVE_SOLVER":
                self.iterative_solver.set_params({item: params[item]})
            elif item.split("-")[1] == "MIXER":
                self.mixer.set_params({item: params[item]})
            elif item.split("-")[1] == "PARAMETERS":
                self.parameters.set_params({item: params[item]})
            else:
                pass
