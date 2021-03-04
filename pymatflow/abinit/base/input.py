# ==============================================================================
# pymatflow.abinit.base.input:
# the general abstract of input file for abinit
# ==============================================================================
import os
import shutil

from pymatflow.abinit.base.electrons import AbinitElectrons
from pymatflow.abinit.base.ions import AbinitIons
from pymatflow.abinit.base.dfpt import AbinitDfpt
from pymatflow.abinit.base.system import AbinitSystem
from pymatflow.abinit.base.properties import AbinitProperties
from pymatflow.abinit.base.guard import AbinitGuard
from pymatflow.abinit.base.misc import AbinitMisc

class AbinitInput:
    """
    """
    def __init__(self):
        self.system = AbinitSystem()
        self.electrons = AbinitElectrons()
        self.ions = AbinitIons()
        self.dfpt = AbinitDfpt()
        self.properties = AbinitProperties()
        self.guard = AbinitGuard()
        self.misc = AbinitMisc()

        self.n = 0 # this is used to control the construction of input string

    def to_string(self):
        """
        :return input_str is the string of all the set params
        """
        # use guard the keep safe
        #self.guard.check_all(system=self.system, electrons=self.electrons, ions=self.ions, dfpt=self.dfpt)
        input_str = ""
        if self.system.status == True:
            input_str += self.system.to_string(n=self.n)
        if self.electrons.status == True:
            input_str += self.electrons.to_string(n=self.n)
        if self.ions.status == True:
            input_str += self.ions.to_string(n=self.n)
        if self.dfpt.status == True:
            input_str += self.dfpt.to_string(n=self.n)
        if self.properties.status == True:
            input_str += self.properties.to_string(n=self.n)

        if self.misc.status == True:
            input_str += self.misc.to_string(n=self.n)
        return input_str

    def get_xyz(self, xyzfile):
        self.system.xyz.get_xyz(xyzfile)

    def set_params(self, params={}):
        for item in params:
            if item in self.electrons.incharge:
                #self.electrons.params[item] = params[item]
                self.electrons.set_param(item, params[item])
                continue
            elif item in self.ions.incharge:
                self.ions.set_param(item, params[item])
                #self.ions.params[item] = params[item]
                continue
            elif item in self.dfpt.incharge:
                #self.dfpt.params[item] = params[item]
                self.dfpt.set_param(item, params[item])
            elif item in self.properties.incharge:
                #self.properties.params[item] = params[item]
                self.properties.set_param(item, params[item])
            else:
                #self.misc.params[item] = params[item]
                self.misc.set_param(item, params[item])

    def set_kpoints(self, kpoints={}):
        #self.electrons.kpoints.set_params(kpoints)
        for item in kpoints:
            if item in self.electrons.kpoints.incharge:
                self.electrons.kpoints.set_param(item, kpoints[item])

    def set_properties(self, properties=[]):
        self.properties.get_option(option=properties)
    #

    def dft_plus_u(self):
        self.electrons.dft_plus_u()
