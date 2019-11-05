#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

import argparse
from emuhelper.cp2k.md import md_run

"""
Usage:
    python aimd_cp2k.py xxx.xyz 
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.
"""

force_eval = {}
motion = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of the calculation", type=str, default="tmp-cp2k-md")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--cutoff", help="CUTOFF, default value: 100 Ry", type=int, default=100)
    parser.add_argument("--xc-functional", help="XC_FUNCTIONAL type", type=str, default="PBE")
    parser.add_argument("--rel-cutoff", help="REL_CUTOFF, default value: 60 Ry", type=int , default=60)
    parser.add_argument("-k", "--kpoints", help="set kpoints like '3 3 3 0 0 0'", type=str, default="3 3 3 0 0 0")
    parser.add_argument("--added-mos", help="ADDED_MOS for SCF", type=int, default=0)
    parser.add_argument("--smear", help="smearing type", type=str, default="FERMI_DIRAC")
    parser.add_argument("--electronic-temp", help="ELECTRON_TEMPERATURE for FERMI_DIRAC SMEAR", type=float, default=300)
    parser.add_argument("--window-size", help="Size of the energy window centred at the Fermi level for ENERGY_WINDOW type smearing", type=float, default=0)

    # motion related parameters
    parser.add_argument("--md-steps", type=int, default=20,
            help="MOTION/MD/STEPS")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]
    force_eval["CUTOFF"] = args.cutoff
    force_eval["XC_FUNCTIONAL"] = args.xc_functional
    force_eval["ADDED_MOS"] = args.added_mos
    force_eval["SMEAR"] = args.smear
    force_eval["ELECTRONIC_TEMPERATURE"] = args.electronic_temp
    force_eval["WINDOW_SIZE"] = args.window_size
    motion["STEPS"] = args.md_steps

    task = md_run(args.file)
    task.md(directory=args.directory, runopt="genrun", force_eval=force_eval, motion=motion)
    task.analysis(directory=args.directory)

