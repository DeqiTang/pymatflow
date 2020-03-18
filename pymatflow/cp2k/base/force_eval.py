"""
a representation for FORCE_EVAL
"""
import numpy as np
import sys
import os
import shutil

from pymatflow.cp2k.base.bsse import cp2k_bsse
from pymatflow.cp2k.base.dft import cp2k_dft
from pymatflow.cp2k.base.eip import cp2k_eip
from pymatflow.cp2k.base.embed import cp2k_embed
from pymatflow.cp2k.base.external_potential import cp2k_external_potential
from pymatflow.cp2k.base.mixed import cp2k_mixed
from pymatflow.cp2k.base.mm import cp2k_mm
from pymatflow.cp2k.base.print import cp2k_print
from pymatflow.cp2k.base.properties import cp2k_properties
from pymatflow.cp2k.base.pw_dft import cp2k_pw_dft
from pymatflow.cp2k.base.qmmm import cp2k_qmmm
from pymatflow.cp2k.base.rescale_forces import cp2k_rescale_forces
from pymatflow.cp2k.base.subsys import cp2k_subsys


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
    def __init__(self):
        self.params = {
                "METHOD": "QS",
                "EMBED": None,
                "STRESS_TENSOR": None,
                }
        self.status = False

        self.subsys = cp2k_subsys()
        self.dft = cp2k_dft()
        self.properties = cp2k_properties()


    def to_input(self, fout):
        # check before write input file
        self.check_spin()

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
            elif item.split("-")[0] == "SUBSYS":
                self.subsys.set_params({item: params[item]})
            elif item.split("-")[0] == "DFT":
                self.dft.set_params({item: params[item]})
            elif item.split("-")[0] == "PROPERTIES":
                self.properties.set_params({item: params[item]})
            else:
                pass
