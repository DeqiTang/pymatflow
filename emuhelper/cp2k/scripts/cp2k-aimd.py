#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from emuhelper.cp2k.md import md_run

"""
Usage:
    python aimd_cp2k.py xxx.xyz 
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.
"""

# =======================================
#           Control Parameters
# =======================================

if __name__ == '__main__':
    task = md_run(sys.argv[1])
    task.ir_spectra()
    task.gen_input()
    task.run() 
    task.analysis()
