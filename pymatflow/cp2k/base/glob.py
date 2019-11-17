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

class cp2k_glob:
    """

    """
    def __init__(self, project_name="ab-initio", run_type="NONE"):
        self.params = {
                "PROJECT": project_name,
                "RUN_TYPE": run_type,
                "PRINT_LEVEL": "LOW",
                "FFTW_PLAN_TYPE": "ESTIMATE",
                }

    def to_input(self, fout):
        fout.write("&GLOBAL\n")
        for item in self.params:
            fout.write("\t%s %s\n" % (item, self.params[item]))
        fout.write("&END GLOBAL\n")
        fout.write("\n")

    def basic_setting(self, run_type="ENERGY_FORCE"):
        """
        """
        if run_type == "ENERGY_FORCE":
            self.params["RUN_TYPE"] = "ENERGY_FORCE"
        elif run_type == "VIBRATIONAL_ANALYSIS":
            self.params["RUN_TYPE"] = "VIBRATIONAL_ANALYSIS"
        elif run_type == "MD":
            self.params["RUN_TYPE"] = "MD"
        elif run_type == "LINEAR_RESPONSE":
            self.params["RUN_TYPE"] = "LINEAR_RESPONSE"
