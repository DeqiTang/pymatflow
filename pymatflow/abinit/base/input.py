
import os
import shutil

from pymatflow.abinit.base.electrons import abinit_electrons
from pymatflow.abinit.base.ions import abinit_ions
from pymatflow.abinit.base.dfpt import abinit_dfpt
from pymatflow.abinit.base.system import abinit_system
from pymatflow.abinit.base.properties import abinit_properties
from pymatflow.abinit.base.guard import abinit_guard


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

        self.electrons.basic_setting()

    def to_input(self, fout):
        self.system.to_input(fout)
        self.electrons.to_input(fout)
        #self.ions.to_input(fout)
        #self.dfpt.to_input(fout)
        self.properties.to_input(fout)


    def get_xyz(self, xyzfile):
        self.system.xyz.get_xyz(xyzfile)

    def set_params(self, electrons={}, ions={}):
        self.electrons.set_params(electrons)
        self.ions.set_params(ions)

    def set_kpoints(self, kpoints={}):
        self.electrons.kpoints.set_params(kpoints)

    def set_properties(self, properties=[]):
        self.properties.get_option(option=properties)

    #

    def dft_plus_u(self):
        self.electrons.dft_plus_u()
