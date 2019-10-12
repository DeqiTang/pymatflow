#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

import pymatgen as mg

"""
usage:
"""

class qe_electrons:
    """

    """
    def __init__(self):
        self.params = {
                "electron_maxstep": None,
                "scf_must_converge": None,
                "conv_thr": None,
                "adaptive_thr": None,
                "conv_thr_init": None,
                "conv_thr_multi": None,
                "mixing_mode": None,
                "mixing_beta": None,
                "mixing_ndim": None,
                "mixing_fixed_ns": None,
                "diagonalization": None,
                "ortho_para": None,
                "diago_thr_init": None,
                "diago_cg_maxiter": None,
                "diago_david_ndim": None,
                "diago_full_acc": None,
                "efield": None,
                "efield_cart": None,
                "efield_phase": None,
                "startingpot": None,
                "startingwfc": None,
                "tqr": None,
                "real_space": None,
                }
    def to_in(self, fout):
        # fout: a file stream for writing
        fout.write("&electrons\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s = %s\n" % (item, str(self.params[item])))
        fout.write("/\n")
        fout.write("\n")
