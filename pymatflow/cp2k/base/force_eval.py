#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from pymatflow.cp2k.base.subsys import cp2k_subsys
from pymatflow.cp2k.base.dft import cp2k_dft
from pymatflow.cp2k.base.properties import cp2k_properties


"""
Usage:
"""

class cp2k_force_eval:
    """
    Note:
        cp2k_force_eval is the class to be responsible for CP2K/FORCE_EVAL
        related settings.

        FORCE_EVAL manage parameters that are needed to describe your system
        and calculate energy and energy of the defined system. so it is actually
        the core driver for other calculations, and as a result, we must set 
        appropriate parameters here to guarantee a proper running for both
        static calculation and other calculations including geometric optimization,
        molecular dynamics, nudged elastic band, etc.

        the relationship between cp2k_force_eval and cp2k_dft is 'have a',
        rather than 'is kind of a', so here combination is used to organize
        the these classes.
    """
    def __init__(self, xyz_f):
        self.params = {
                "METHOD": "QS",
                "EMBED": None,
                "STRESS_TENSOR": None,
                }
        self.status = False

        self.subsys = cp2k_subsys(xyz_f)
        self.dft = cp2k_dft()
        self.properties = cp2k_properties()

        self.check_spin()

    def to_input(self, fout):
        fout.write("&FORCE_EVAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t%s %s\n" % (item, self.params[item]))
        if self.subsys.status == True:
            self.subsys.to_input(fout)
        if self.dft.status == True:
            self.dft.to_input(fout)
        if self.properties.status == True:
            self.properties.to_input(fout)
        fout.write("&END FORCE_EVAL\n") 
        fout.write("\n")
    
    def check_spin(self):
        """
        call self.dft.check_spin()
        """
        self.dft.check_spin(self.subsys.xyz)

    def basic_setting(self):
        self.subsys.status = True
        self.dft.status = True
        self.dft.mgrid.params["CUTOFF"] = 100
        self.dft.mgrid.params["REL_CUTOFF"]= 60

    def set_params(self, params):
        """
        parameters for sub section(like dft), are handled over 
        to sub section controllers.
        """
        for item in params:
            if len(item.split("-")) == 1:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[0] == "DFT":
                self.dft.set_params({item: params[item]})
