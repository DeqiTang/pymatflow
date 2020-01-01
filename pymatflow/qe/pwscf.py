#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.qe.base.control import qe_control
from pymatflow.qe.base.system import qe_system
from pymatflow.qe.base.electrons import qe_electrons
from pymatflow.qe.base.ions import qe_ions
from pymatflow.qe.base.cell import qe_cell
from pymatflow.qe.base.arts import qe_arts


class pwscf:
    """
    About:
        this class provide a base representation of pwscf
    """
    def __init__(self):
        self.control = qe_control()
        self.system = qe_system()
        self.electrons = qe_electrons()
        self.ions = qe_ions()
        self.cell = qe_cell()
        self.arts = qe_arts()
        
        #self.control.basic_setting("scf") 
        self.electrons.basic_setting()
        self.set_kpoints() # default kpoint setting

    def get_xyz(self, xyzfile):
        """
        xyz_f:
            a modified xyz formatted file(the second line specifies the cell of the 
            system).
        """
        self.arts.xyz.get_xyz(xyzfile)
        self.system.basic_setting(self.arts)
        self.arts.basic_setting(ifstatic=True)


    def set_params(self, control={}, system={}, electrons={}, ions={}, cell={}):
        # check if user try to set occupations and smearing and degauss
        # through system. if so, use self.set_occupations() which uses
        # self.system.set_occupations() to set them, as self.system.set_params() 
        # is suppressed from setting occupations related parameters
        self.set_occupations(system)
        self.control.set_params(control)
        self.system.set_params(system)
        self.electrons.set_params(electrons)
        self.ions.set_params(ions)
        self.cell.set_params(cell)

    def set_kpoints(self, kpoints_option="automatic", kpoints_mp=[2, 2, 2, 0, 0, 0]):
        self.arts.set_kpoints(option=kpoints_option, kpoints_mp=kpoints_mp)
 
    def set_atomic_forces(self, pressure=None, pressuredir=None):
        self.arts.set_atomic_forces(pressure=pressure, direction=pressuredir)       

    def set_occupations(self, system):
        """
            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.system.set_occupations() to 
            # set them, as self.system.set_params() is suppressed from setting
            # occupations related parameters
            # if occupations == None, use default smearing occupation. and 
            # if occupations == "tetrahedra" the value set for smearing and degauss is ignored.
            # if occupations == "smearing", the value of smearing and degauss
            # should be legacy, not None or other illegal values.
        """
        if "occupations" in system:
            if system["occupations"] == None: # user default setting of set_occupations()
                self.system.set_occupations()
            elif system["occupations"] == "tetrahedra":
                self.system.set_occupations(occupations="tetrahedra")
            elif system["occupations"] == "smearing":
                if "smearing" in system and "degauss" in system:
                    self.system.set_occupations(occupations="smearing", smearing=system["smearing"], degauss=system["degauss"])
                elif "smearing" in system:
                    self.system.set_occupations(occupations="smearing", smearing=system["smearing"])
                elif "degauss" in system:
                    self.system.set_occupations(occupations="smearing", degauss=system["degauss"])
                else:
                    self.system.set_occupations(occupations="smearing")
            elif system["occupations"] == "tetrahedra_lin":
                self.system.set_occupations(occupations="tetrahedra_lin")
            elif system["occupations"] == "tetrahedra_opt":
                self.system.set_occupations(occupations="tetrahedra_opt")
            elif system["occupations"] == "fixed":
                self.system.set_occupations(occupations="fixed")
            elif system["occupations"] == "from_input":
                self.system.set_occupations(occupations="from_input")

    def gen_yh(self, inpname, output, directory, cmd="pw.x"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".in")[0]+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % (cmd, inpname, output))

