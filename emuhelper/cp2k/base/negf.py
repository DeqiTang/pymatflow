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

class cp2k_negf:
    """

    """
    def __init__(self):
        self.params = {
                "DELTA_NPOLES": None,
                "DISABLE_CACHE": None,
                "ENERGY_LBOUND": None,
                "EPS_DENSITY": None,
                "EPS_GEO": None,
                "EPS_GREEN": None,
                "EPS_SCF": None,
                "ETA": None,
                "GAMMA_KT": None,
                "HOMO_LUMO_GAP": None,
                "INTEGRATION_MAX_POINTS": None,
                "INTEGRATION_METHOD": None,
                "INTEGRATION_MIN_POINTS": None,
                "MAX_SCF": None,
                "NPROC_POINT": None,
                "V_SHIFT": None,
                "V_SHIFT_MAX_ITERS": None,
                "V_SHIFT_OFFSET": None,
                }

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("&NEGF\n")
        for item in self.params:
            fout.write("\t%s %s\n" % (item, self.params[item]))
        fout.write("&END NEGF\n")
        fout.write("\n")
