#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

"""
Usage:
"""
# ==================================
# CP2K / FORCE_EVAL / DFT / ALMO_SCF
# ==================================
class cp2k_dft_almo_scf:
    def __init__(self):
        pass

# ============================================
# ============================================
class cp2k_dft_auxiliary_density_matrix_method:
    def __init__(self):
        pass

# ================================
# ================================
class cp2k_dft_density_fitting:
    def __init__(self):
        pass

# ====================
# ====================
class cp2k_dft_efield:
    def __init__(self):
        self.params = {
                "ENVELOP": None,
                "INTENSITY": None,
                "PHASE": None,
                "POLARISATION": None,
                "WAVELENGTH": None,
                }

# ==============================
# ==============================
class cp2k_dft_external_density:
    def __init__(self):
        pass

# ==================================
# ==================================
class cp2k_dft_external_potential:
    def __init__(self):
        pass

# ==========================
# ==========================
class cp2k_dft_external_vxc:
    def __init__(self):
        pass

# ==========================
# ==========================
class cp2k_dft_kg_method:
    def __init__(self):
        pass

# =================================
# CP2K / FORCE_EVAL / DFT / KPOINTS
# =================================
class cp2k_dft_kpoints:
    def __init__(self):
        self.params = {
                "EPS_GEO": None,
                "FULL_GRID": None,
                "KPOINT": None,
                "PARALLEL_GROUP_SIZE": None,
                "SCHEME": None,
                "SYMMETRY": None,
                "UNITS": None,
                "VERBOSE": None,
                "WAVEFUNCTIONS": None,
                }
    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&KPOINTS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END KPOINTS\n")

# =================================
# CP2K / FORCE_EVAL /DFT / LOCALIZE
# =================================
class cp2k_dft_localize:
    """
    About:
        to calculate IR spectra from md running
        it is necessary to have dipole information
        for the molecules available in the simulated
        trajectory.
    """
    def __init__(self):
        self.params = {
                "METHOD": None,
                "MAX_ITER": None,
                }
        self.status = False

    def to_dft(self, fout):
        # fout: a file stream for writing
        if self.status == False:
            return
        fout.write("\t\t&LOCALIZE %s\n" % ".TRUE.")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&PRINT\n")
        fout.write("\t\t\t\t&WANNIER_CENTERS\n")
        fout.write("\t\t\t\t\tIONS+CENTERS\n")
        fout.write("\t\t\t\t\tFILENAME xxx_wannier.xyz\n")
        fout.write("\t\t\t\t&END WANNIER_CENTERS\n")
        fout.write("\t\t\t&END PRINT\n")
        fout.write("\t\t&END LOCALIZE\n")


# =======================================
# CP2K / FORCE_EVAL / DFT / LOW_SPIN_ROKS
# =======================================
class cp2k_dft_low_spin_roks:
    def __init__(self):
        self.params = {
                "ENERGY_SCALING": None,
                "SPIN_CONFIGURATION": None
                }
    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&LOW_SPIN_ROKS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END LOW_SPIN_ROKS\n")

# =======================================
# CP2K / FORCE_EVAL / DFT /  LS_SCF
# =======================================
class cp2k_dft_ls_scf:
    def __init__(self):
        self.params = {
                "CHECK_S_INV": None,
                "DYNAMIC_THRESHOLD": None,
                "EPS_DIIS": None,
                "EPS_FILTER": None,
                "EPS_LANCZOS": None,
                "EPS_SCF": None,
                "EXTRAPOLATION_ORDER": None,
                "FIXED_MU": None,
                "INI_DIIS": None,
                "LS_DIIS": None,
                "MATRIX_CLUSTER_TYPE": None,
                "MAX_DIIS": None,
                "MAX_ITER_LANCZOS": None,
                "MAX_SCF": None,
                "MIXING_FRACTION": None,
                "MU": None,
                "NMIXING": None,
                "NON_MONOTONIC": None,
                "PERFORM_MU_SCAN": None,
                "PURIFICATION_METHOD": None,
                "REPORT_ALL_SPARSITIES": None,
                "RESTART_READ": None,
                "RESTART_WRITE": None,
                "SIGN_SQRT_ORDER": None,
                "SINGLE_PRECISION_MATRICES": None,
                "S_INVERSION": None,
                "S_PRECONDITIONER": None,
                }
    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&LS_SCF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END LS_SCF\n")
 
# ===============================
# CP2K / FORCE_EVAL / DFT / MGRID
# ===============================
class cp2k_dft_mgrid:
    def __init__(self):
        self.params = {
                "CUTOFF": 280,
                "REL_CUTOFF": 40,
                "NGRIDS": 4,
                "COMMENSURATE": None,
                "MULTIGRID_CUTOFF": None,
                "MULTIGRID_SET": None,
                "PROGRESSION_FACTOR": None,
                "SKIP_LOAD_BALANCE_DISTRIBUTED": None,
                }
    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&MGRID\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END MGRID\n")


# ================================
# ================================
class cp2k_dft_periodic_efield:
    def __init__(self):
        pass

# =================================
# CP2K / FORCE_EVAL / DFT / POISSON
# =================================
class cp2k_dft_poisson:
    def __init__(self):
        self.params = {
                "PERIODIC": "XYZ",
                "POISSON_SOLVER": "PERIODIC",
                }
    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&POISSON\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END POISSON\n")

# ===============================
# CP2K / FORCE_EVAL / DFT / PRINT
# ===============================
class cp2k_dft_print:
    def __init__(self):
        self.params = {
                }
        self.electron_density = False
        self.bands = False
        self.pdos = False
    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.electron_density == True:
            fout.write("\t\t\t&E_DENSITY_CUBE\n")
            fout.write("\t\t\t&END E_DENSITY_CUBE\n")
        if self.bands == True:
            fout.write("\t\t\t&BAND_STRUCTURE\n")
            fout.write("\t\t\t\tADDED_MOS 10\n") 
            fout.write("\t\t\t\tFILE_NAME bands.bs\n")
            fout.write("\t\t\t\t&KPOINT_SET\n")
            fout.write("\t\t\t\t\tUNITS B_VECTOR\n")
            fout.write("\t\t\t\t\tSPECIAL_POINT GAMA 0.0000000000     0.0000000000     0.0000000000\n") # GAMA
            fout.write("\t\t\t\t\tSPECIAL_POINT X    0.5000000000     0.0000000000     0.5000000000\n") # X
            fout.write("\t\t\t\t\tSPECIAL_POINT W    0.5000000000     0.2500000000     0.7500000000\n") # W
            fout.write("\t\t\t\t\tSPECIAL_POINT GAMA 0.0000000000     0.0000000000     0.0000000000\n") # GAMA
            fout.write("\t\t\t\t\tSPECIAL_POINT U    0.6250000000     0.2500000000     0.6250000000\n") # U
            fout.write("\t\t\t\t\tSPECIAL_POINT W    0.5000000000     0.2500000000     0.7500000000\n") # W
            fout.write("\t\t\t\t\tNPOINTS 20\n")
            fout.write("\t\t\t\t&END KPOINT_SET\n")
            fout.write("\t\t\t&END BAND_STRUCTURE\n")
        if self.pdos == True:
            fout.write("\t\t\t&PDOS\n")
            fout.write("\t\t\t\tNLUMO -1\n") # print all projected DOS available
            fout.write("\t\t\t\tCOMPONENTS\n") # split the density by quantum number
            fout.write("\t\t\t&END PDOS\n")
        fout.write("\t\t&END PRINT\n")
    def print_electron_density(self):
        self.electron_density = True

    def print_bands(self):
        self.bands = True
    def print_pdos(self):
        self.pdos = True

# ============================
# CP2K / FORCE_EVAL / DFT / QS
# ============================
class cp2k_dft_qs:
    def __init__(self):
        self.params = {
                "METHOD": "GPW",
                "EPS_DEFAULT": "1.0E-10",
                "FORCE_PAW": None,
                }
    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&QS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END QS\n")


# ==================================
# ==================================
class cp2k_dft_real_time_propagation:
    def __init__(self):
        self.params = {
                "ACCURACY_REFINEMENT": None,
                "APPLY_DELTA_PULSE": None,
                "ASPC_ORDER": None,
                "PERIODIC": None,
                "PROPAGATOR": None,
                "MAX_EXP": None,
                "MAX_ITER": None,
                "EPS_ITER": None,
                }


# ==========================
# ==========================
class cp2k_dft_relativistic:
    def __init__(self):
        self.params = {
                "DKH_ORDER": None,
                "METHOD": None,
                "POTENTIAL": None,
                "TRANSFORMATION": None,
                "ZORA_TYPE": None,
                "Z_CUTOFF": None,
                }

# ==================
# ==================
class cp2k_dft_sccs:
    def __init__(self):
        pass

# =============================
# CP2K / FORCE_EVAL / DFT / SCF
# =============================
class cp2k_dft_scf_diagonalization:
    def __init__(self):
        self.section = '.TRUE.'
        self.params = {
                "ALGORITHM": "STANDARD",
                "EPS_ADAPT": None,
                "EPS_ITER": None,
                "EPS_JACOBI": None,
                "JACOBI_THRESHOLD": None,
                "MAX_ITER": None,
                }
    def to_scf(self, fout):
       """
       fout: a file stream for writing
       """
       fout.write("\t\t\t&DIAGONALIZATION\n")
       for item in self.params:
           if self.params[item] is not None:
               fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
       fout.write("\t\t\t&END DIAGONALIZATION\n")

class cp2k_dft_scf_mixing:
    def __init__(self):
        self.section = '.TRUE.'
        self.params = {
                "ALPHA": 0.4,
                "BETA": None,
                "METHOD": "BROYDEN_MIXING",
                }
    def to_scf(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&MIXING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END MIXING\n")

class cp2k_dft_scf_smear:
    def __init__(self):
        self.section = '.TRUE.'
        self.params = {
                "METHOD": 'FERMI_DIRAC',
                "ELECTRONIC_TEMPERATURE": None,
                }
    def to_scf(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&SMEAR %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END SMEAR\n")

class cp2k_dft_scf_print:
    def __init__(self):
        self.params = {
                "DM_RESTART_WRITE": None,
                }
    def to_scf(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END PRINT\n")


class cp2k_dft_scf:
    def __init__(self):
        self.params = {
                "ADDED_MOS": 20,
                "SCF_GUESS": "ATOMIC",
                "EPS_SCF": "1.0E-06",
                "MAX_SCF": 50,
                "MAX_DIIS": None,
                "MAX_ITER_LUMO": None,
                "MAX_SCF_HISTORY": None,
                "ROKS_SCHEME": None,
                }
        self.diagonalization = cp2k_dft_scf_diagonalization()
        self.mixing = cp2k_dft_scf_mixing()
        self.smear = cp2k_dft_scf_smear()
        self.printout = cp2k_dft_scf_print()
        self.ifsmear = True

    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&SCF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        self.diagonalization.to_scf(fout)
        self.mixing.to_scf(fout)
        if self.ifsmear == True:
            if self.params["ADDED_MOS"] == None or self.params["ADDED_MOS"] == 0:
                print("If you are using smearing, you should set ADDED_MOS too!!!\n")
                sys.exit()
            self.smear.to_scf(fout)
        self.printout.to_scf(fout)
        fout.write("\t\t&END SCF\n")


# ===================
# ===================
class cp2k_dft_scrf:
    def __init__(self):
        pass


# ===================
# ===================
class cp2k_dft_sic:
    def __init__(self):
        self.params = {
                "ORBITAL_SET": None,
                "SIC_METHOD": None,
                "SIC_SCALING_A": None,
                "SIC_SCALING_B": None,
                }       


# ====================
# ====================
class cp2k_dft_tddfpt:
    def __init__(self):
        self.params = {
                "CONVERGENCE": None,
                "DIAG_METHOD": None,
                "INVERT_S": None,
                "KERNEL": None,
                "LSD_SINGLETS": None,
                "MAX_KV": None,
                "NEV": None,
                "NLUMO": None,
                "NREORTHO": None,
                "OE_CORR": None,
                "PRECONDITIONER": None,
                "RESTARTS": None,
                "RES_ETYPE": None,
                }       

# ==============================
# ==============================
class cp2k_dft_transport:
    def __init__(self):
        self.params = {
                "TRANSPORT_METHOD": None,
                }       
# ==============================
# ==============================
class cp2k_dft_xas():
    def __init__(self):
        self.params = {
                "ADDED_MOS": None,
                "METHOD": None,
                }

# =============================
# CP2K / FORCE_EVAL / DFT / XC
# =============================
class cp2k_dft_xc_xc_functional:
    def __init__(self):
        self.section = "PBE"
        self.params = {
                }
    def to_xc(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&XC_FUNCTIONAL %s\n" % self.section)
        fout.write("\t\t\t&END XC_FUNCTIONAL\n")


class cp2k_dft_xc:
    def __init__(self):
        self.params = {
                
                }
        self.xc_functional = cp2k_dft_xc_xc_functional()
    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&XC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s" % (item, str(self.params[item])))
        self.xc_functional.to_xc(fout)
        fout.write("\t\t&END XC\n")


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
        self.mgrid.to_dft(fout)
        self.xc.to_dft(fout)
        #self.kpoints.to_dft(fout)
        self.scf.to_dft(fout)
        self.localize.to_dft(fout)
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
