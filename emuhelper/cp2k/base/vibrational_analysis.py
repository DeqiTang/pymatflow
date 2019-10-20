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
    Note:
        if you are calculating IR, you should print out
        dipole moment throught DFT/PRINT/MOMENTS.
        for direct diagonalization method dipole moment
        calculation of Berry phase approach is not supported
        , and for OT method, we can use PERIODIC .TRUE in
        DFT/PRINT/MOMENTS

        inf DFT/PRINT/MOMENTS
        PERIODIC:
            Use Berry phase formula (PERIODIC=T) or simple operator (PERIODIC=F). 
            The latter normally requires that the CELL is periodic NONE.
        so if we want to use Berry phase to calculate dipole moment 
        we have to use OT in SCF, and if we use PERIODIC=F, we have to set 
        PERIODIC=NONE in CELL.
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
            if self.params[item] is not None:
                fout.write("\t%s %s\n" % (item, self.params[item]))
        fout.write("&END VIBRATIONAL_ANALYSIS\n")
        fout.write("\n")

    def basic_setting(self):
        self.params["INTENSITIES"] = True
        self.params["DX"] = 1.0e-2 # default value in Bohr
        self.params["TC_PRESSURE"] = 1.01325e5  # Pa
        self.params["TC_TEMPERATURE"] = 2.7315e2 # K
        self.params["THERMOCHEMISTRY"] = False # Calculation of the thermochemical data. Valid for molecules in the gas phase.
