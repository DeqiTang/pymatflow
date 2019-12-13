#!/usr/bin/env python
# _*_ coding: utf-8 _*_

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
        self.status = False

        self.params["EPS_FILTER"] = 1.0E-7
        self.params["EPS_SCF"] = 1.0E-5
        self.params["S_PRECONDITIONER"] = "ATOMIC"

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&LS_SCF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END LS_SCF\n")

