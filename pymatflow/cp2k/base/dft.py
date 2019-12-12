#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from pymatflow.cp2k.base.dft_almo_scf import cp2k_dft_almo_scf
from pymatflow.cp2k.base.dft_auxiliary_density_matrix_method import cp2k_dft_auxiliary_density_matrix_method
from pymatflow.cp2k.base.dft_density_fitting import cp2k_dft_density_fitting
from pymatflow.cp2k.base.dft_external_density import cp2k_dft_external_density
from pymatflow.cp2k.base.dft_external_potential import cp2k_dft_external_potential
from pymatflow.cp2k.base.dft_external_vxc import cp2k_dft_external_vxc
from pymatflow.cp2k.base.dft_kg_method import cp2k_dft_kg_method
from pymatflow.cp2k.base.dft_low_spin_roks import cp2k_dft_low_spin_roks
from pymatflow.cp2k.base.dft_ls_scf import cp2k_dft_ls_scf
from pymatflow.cp2k.base.dft_mgrid import cp2k_dft_mgrid
from pymatflow.cp2k.base.dft_periodic_efield import cp2k_dft_periodic_efield
from pymatflow.cp2k.base.dft_poisson import cp2k_dft_poisson
from pymatflow.cp2k.base.dft_qs import cp2k_dft_qs
from pymatflow.cp2k.base.dft_real_time_propagation import cp2k_dft_real_time_propagation
from pymatflow.cp2k.base.dft_relativistic import cp2k_dft_relativistic
from pymatflow.cp2k.base.dft_sccs import cp2k_dft_sccs
from pymatflow.cp2k.base.dft_scf import cp2k_dft_scf
from pymatflow.cp2k.base.dft_scrf import cp2k_dft_scrf
from pymatflow.cp2k.base.dft_sic import cp2k_dft_sic
from pymatflow.cp2k.base.dft_tddfpt import cp2k_dft_tddfpt
from pymatflow.cp2k.base.dft_transport import cp2k_dft_transport
from pymatflow.cp2k.base.dft_xas import cp2k_dft_xas
from pymatflow.cp2k.base.dft_xc import cp2k_dft_xc
from pymatflow.cp2k.base.dft_print import cp2k_dft_print
from pymatflow.cp2k.base.dft_kpoints import cp2k_dft_kpoints
from pymatflow.cp2k.base.dft_localize import cp2k_dft_localize
from pymatflow.cp2k.base.dft_efield import cp2k_dft_efield

"""
usage:
"""

# ============================================
# CP2K / DFT
#=============================================

class cp2k_dft:
    """

    """
    def __init__(self):
        """
        BASIS_MOLOPT含有所有元素的DZVP-MOLOPT-SR-GTH基组
        GTH_POTENTIALS含有所有元素的GTH-PBE赝势以及几乎所有
        元素的GTH-BLYP赝势(似乎Nb没有). 所以将其设为默认值
        """
        self.params = {
                "AUTO_BASIS": None,
                "BASIS_SET_FILE_NAME": "BASIS_MOLOPT",
                "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
                "CHARGE": None,
                "EXCITATIONS": None,
                "MULTIPLICITY": None,
                "PLUS_U_METHOD": None,
                "RELAX_MULTIPLICITY": None,
                "ROKS": None,
                "SUBCELLS": None,
                "SURFACE_DIPOLE_CORRECTION": None,
                "SURF_DIP_DIR": None,
                "LSD": None, # alis: LSD = SPIN_POLARIZED = UNRESTRICTED_KOHN_SHAM = UKS
                "WFN_RESTART_FILE_NAME": None,
                }

        self.almo_scf = cp2k_dft_almo_scf()

        self.auxiliary_density_matrix_method = cp2k_dft_auxiliary_density_matrix_method()

        self.density_fitting = cp2k_dft_density_fitting()

        self.efield = cp2k_dft_efield() 
        
        self.external_density = cp2k_dft_external_density()

        self.external_potential = cp2k_dft_external_potential()

        self.external_vxc = cp2k_dft_external_vxc()

        self.kg_method = cp2k_dft_kg_method()

        self.kpoints = cp2k_dft_kpoints()

        self.localize = cp2k_dft_localize()

        self.low_spin_roks = cp2k_dft_low_spin_roks()

        self.ls_scf = cp2k_dft_ls_scf()
        
        self.mgrid = cp2k_dft_mgrid()

        self.periodic_efield = cp2k_dft_periodic_efield()

        self.poisson = cp2k_dft_poisson()
        
        self.printout = cp2k_dft_print()

        self.qs = cp2k_dft_qs()
        
        self.real_time_propagation = cp2k_dft_real_time_propagation()
        
        self.relativistic = cp2k_dft_relativistic() 
        
        self.sccs = cp2k_dft_sccs()

        self.scf = cp2k_dft_scf()
        
        self.scrf = cp2k_dft_scrf()
        
        self.sic = cp2k_dft_sic() 
        
        self.tddfpt = cp2k_dft_tddfpt()
        
        self.transport = cp2k_dft_transport()
        
        self.xas = cp2k_dft_xas() 
        
        self.xc = cp2k_dft_xc()

    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&DFT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        self.qs.to_dft(fout)
        self.poisson.to_dft(fout)
        if self.ls_scf.section.upper() == "TRUE":
            self.ls_scf.to_dft(fout)
        self.mgrid.to_dft(fout)
        self.xc.to_dft(fout)
        self.kpoints.to_dft(fout)
        self.scf.to_dft(fout)
        self.localize.to_dft(fout)
        if self.periodic_efield.section.upper() == "TRUE":
            self.periodic_efield.to_dft(fout)
        self.printout.to_dft(fout)
        fout.write("\t&END DFT\n")
        #fout.write("\n")

    def check_spin(self, xyz):
        """
        xyz: base_xyz or cp2k_xxyz or cp2k_subsys
        """
        n_electrons = 0
        for atom in xyz.atoms:
            n_electrons += mg.Element(atom.name).number
        if n_electrons % 2 == 1:
            self.params["LSD"] = ".TRUE."
        else:
            self.params["LSD"] = None

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                if item.split("-")[-1] == "LS_SCF":
                    self.ls_scf.section = params[item]
                else:
                    self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "QS":
                self.qs.set_params({item: params[item]})
            elif item.split("-")[1] == "MGRID":
                self.mgrid.set_params({item: params[item]})
            elif item.split("-")[1] == "XC":
                self.xc.set_params({item: params[item]})
            elif item.split("-")[1] == "SCF":
                self.scf.set_params({item: params[item]})
            elif item.split("-")[1] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[1] == "LS_SCF":
                self.ls_scf.set_params({item: params[item]})
            elif item.split("-")[1] == "KPOINTS":
                self.kpoints.set_params({item: params[item]})
