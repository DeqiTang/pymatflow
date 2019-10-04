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

class cp2k_vibrational_analysis:
    """

    """
    def __init__(self):
        self.params = {
                "DX": None,
                "FULLY_PERIODIC": None,
                "INTENSITES": None,
                "NPROC_REP": None,
                "PROC_DIST_TYPE": None,
                "TC_PRESSURE": None,
                "TC_TEMPERATURE": None,
                "THERMOCHEMISTRY": None,
                }

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("&VIBRATIONAL_ANALYSIS\n")
        for item in self.params:
            fout.write("\t%s %s\n" % (item, self.params[item]))
        fout.write("&END VIBRATIONAL_ANALYSIS\n")
        fout.write("\n")

