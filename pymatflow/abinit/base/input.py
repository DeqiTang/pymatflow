# ==============================================================================
# pymatflow.abinit.base.input:
# the general abstract of input file for abinit
# ==============================================================================
import os
import shutil

from pymatflow.abinit.base.electrons import abinit_electrons
from pymatflow.abinit.base.ions import abinit_ions
from pymatflow.abinit.base.dfpt import abinit_dfpt
from pymatflow.abinit.base.system import abinit_system
from pymatflow.abinit.base.properties import abinit_properties
from pymatflow.abinit.base.guard import abinit_guard
from pymatflow.abinit.base.misc import abinit_misc

class abinit_input:
    """
    """
    def __init__(self):
        self.system = abinit_system()
        self.electrons = abinit_electrons()
        self.ions = abinit_ions()
        self.dfpt = abinit_dfpt()
        self.properties = abinit_properties()
        self.guard = abinit_guard()
        self.misc = abinit_misc()

    def to_input(self, fout):
        # use guard the keep safe
        self.guard.check_all(system=self.system, electrons=self.electrons, ions=self.ions, dfpt=self.dfpt)
        #
        self.system.to_input(fout)
        self.electrons.to_input(fout)
        self.ions.to_input(fout)
        self.dfpt.to_input(fout)
        self.properties.to_input(fout)
        self.misc.to_input(fout)

    def to_string(self):
        """
        :return input_str is the string of all the set params
        """
        input_str = ""
        input_str += self.system.to_string()
        input_str += self.electrons.to_string()
        input_str += self.ions.to_string()
        input_str += self.dfpt.to_string()
        input_str += self.properties.to_string()
        input_str += self.misc.to_string()
        return input_str

    def get_xyz(self, xyzfile):
        self.system.xyz.get_xyz(xyzfile)

    def set_params(self, params={}):
        for item in params:
            if item in self.electrons.incharge:
                self.electrons.params[item] = params[item]
                continue
            elif item in self.ions.incharge:
                self.ions.params[item] = params[item]
                continue
            elif item in self.dfpt.incharge:
                self.dfpt.params[item] = params[item]
            elif item in self.properties.incharge:
                self.properties.params[item] = params[item]
            else:
                self.misc.params[item] = params[item]

    def set_kpoints(self, kpoints={}):
        self.electrons.kpoints.set_params(kpoints)

    def set_properties(self, properties=[]):
        self.properties.get_option(option=properties)
    #

    def dft_plus_u(self):
        self.electrons.dft_plus_u()
