#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

import pymatgen as mg

"""
usage:
"""

class qe_system:
    """

    """
    def __init__(self):
        self.params = {
                "ibrav": None,
                "celldm": None,
                "A": None,
                "B": None,
                "C": None,
                "cosAB": None,
                "cosAC": None,
                "cosBC": None,
                "nat": None,
                "ntyp": None,
                "nbnd": None,
                "tot_charge": None,
                "starting_charge": None,
                "tot_magnetization": None,
                "starting_magnetization": None,
                "ecutwfc": None,
                "ecutrho": None,
                "ecutfock": None,
                "nr1": None,
                "nr2": None,
                "nr3": None,
                "nr1s": None,
                "nr2s": None,
                "nr3s": None,
                "nosym": None,
                "nosym_evc": None,
                "noinv": None,
                "no_t_rev": None,
                "force_symmorphic": None,
                "use_all_frac": None,
                "occupations": None,
                "one_atom_occupations": None,
                "starting_spin_angle": None,
                "degauss": None,
                "smearing": None,
                "nspin": None,
                "noncolin": None,
                "ecfixed": None,
                "qcutz": None,
                "q2sigma": None,
                "input_dft": None,
                "exx_fraction": None,
                "screening_parameter": None,
                "exxdiv_treatment": None,
                "x_gamma_extrapolation": None,
                "ecutvcut": None,
                "nqx1": None,
                "nqx2": None,
                "nqx3": None,
                "localization_thr": None,
                "lda_plus_u": None,
                "lda_plus_u_kind": None,
                "Hubbard_U": None,
                "Hubbard_J0": None,
                "Hubbard_alpha": None,
                "Hubbard_beta": None,
                "Hubbard_J": None,
                "starting_ns_eigenvalue": None,
                "U_projection_type": None,
                "edir": None,
                "emaxpos": None,
                "eopreg": None,
                "eamp": None,
                "angle1": None,
                "angle2": None,
                "lforcet": None,
                "constrained_magnetization": None,
                "fixed_magnetization": None,
                "lambda": None,
                "report": None,
                "lspinorb": None,
                "assume_isolated": None,
                "esm_bc": None,
                "esm_w": None,
                "esm_efield": None,
                "esm_nfit": None,
                "fcp_mu": None,
                "vdw_corr": None,
                "london": None,
                "london_s6": None,
                "london_c6": None,
                "london_rvdw": None,
                "london_rcut": None,
                "dftd3_version": None,
                "dftd3_threebody": None,
                "ts_vdw_econv_thr": None,
                "ts_vdw_isolated": None,
                "xdm": None,
                "xdm_a1": None,
                "xdm_a2": None,
                "space_group": None,
                "uniqueb": None,
                "origin_choice": None,
                "rhombohedral": None,
                "zgate": None,
                "relaxz": None,
                "block": None,
                "block_1": None,
                "block_2": None,
                "block_height": None,
                }
    def to_in(self, fout):
        # fout: a file stream for writing
        # ==============================
        # checking legacy of parameters
        # if there is problem with it 
        # the it will kill the script
        # ==============================
        self.check_all()
        # ==============================
        fout.write("&system\n")
        for item in self.params:
            if self.params[item] is not None:
                if type(self.params[item]) == str:
                    fout.write("%s = '%s'\n" % (item, str(self.params[item])))
                else:
                    fout.write("%s = %s\n" % (item, str(self.params[item])))
        fout.write("/\n")
        fout.write("\n")

    def check_all(self):
        """
        """
        must_define = ["ibrav", "nat", "ntyp", "ecutwfc"]
        for item in must_define:
            if self.params[item] is None:
                print("===================================\n")
                print("          Warning !!!!!!\n")
                print("===================================\n")
                print("the following parameters must always be\n")
                print("set:\n")
                print("[%s, %s, %s, %s]\n" % (must_define[0], must_define[1], must_define[2], must_define[3]))
                sys.exit(1)

    def basic_setting(self, arts):
        """
        arts: an object of qe.base.arts.qe_arts
        """
        self.params["ibrav"] = 0
        self.params["nat"] = arts.xyz.natom
        self.params["ntyp"] = arts.xyz.nspecies

        self.params["ecutwfc"] = 100
        self.params["input_dft"] = 'PBE'
        self.params["occupations"] = "smearing"
        self.params["smearing"] = "gaussian"
        self.params["degauss"] = 0.0001
    
    def set_params(self, params):
        """
        params: a dict storing the parameters and values
        """
        for item in params:
            self.params[item] = params[item]
