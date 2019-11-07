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


    parser.add_argument("--ls-scf", type=str, default="FALSE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="use linear scaling scf method")
    
    parser.add_argument("--qs-method", type=str, default="GPW",
            choices=["AM1", "DFTB", "GAPW", "GAPW_XC", "GPW", "LRIGPW", "MNDO", "MNDOD", 
                "OFGPW", "PDG", "PM3", "PM6", "PM6-FM", "PNNL", "RIGPW", "RM1"],
            help="specify the electronic structure method")

    parser.add_argument("--eps-scf", type=float, default=1.0e-5,
            help="DFT-SCF-EPS_SCF")

    parser.add_argument("--xc-functional", type=str, default="PBE",
            help="XC_FUNCTIONAL type")

    parser.add_argument("--cutoff", type=int, default=100,
            help="CUTOFF, default value: 100 Ry")

    parser.add_argument("--rel-cutoff", type=int, default=60,
            help="REL_CUTOFF, default value: 60 Ry")

    parser.add_argument("-k", "--kpoints", help="set kpoints like '3 3 3 0 0 0'", type=str, default="3 3 3 0 0 0")
    
    parser.add_argument("--diag", type=str, default="TRUE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing tranditional diagonalization for SCF")

    parser.add_argument("--ot", type=str, default="FALSE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing orbital transformation for SCF")

    parser.add_argument("--smear", type=str, default="FALSE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="switch on or off smearing for occupation")

    parser.add_argument("--added-mos", help="ADDED_MOS for SCF", type=int, default=0)

    parser.add_argument("--smear-method", help="smearing type: FERMI_DIRAC, ENERGY_WINDOW", type=str, default="FERMI_DIRAC")

    parser.add_argument("--electronic-temp", help="ELECTRON_TEMPERATURE for FERMI_DIRAC SMEAR", type=float, default=300)

    parser.add_argument("--window-size", help="Size of the energy window centred at the Fermi level for ENERGY_WINDOW type smearing", type=float, default=0)


    # motion related parameters
    parser.add_argument("--md-steps", type=int, default=20,
            help="MOTION/MD/STEPS")
    parser.add_argument("--ensemble", type=str, default="NVE",
            choices=["NVE",  "NVT","HYDROSTATICSHOCK", "ISOKIN", "LANGEVIN", "MSST", "MSST_DAMPED"],
            help="MOTION/MD/ENSEMBLE")
    parser.add_argument("--temperature", type=str, default=300,
            help="The temperature in K used to initialize the velocities with init and pos restart, and in the NPT/NVT simulations")
    parser.add_argument("--traj-format", type=str, default="XMOL",
            help="type of output trajectory for MOTION, note: DCD format can be visualized by vmd")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]

    force_eval["DFT-LS_SCF"] = args.ls_scf
    force_eval["DFT-QS-METHOD"] = args.qs_method
    force_eval["DFT-MGRID-CUTOFF"] = args.cutoff
    force_eval["DFT-MGRID-REL_CUTOFF"] = args.rel_cutoff
    force_eval["DFT-XC-XC_FUNCTIONAL"] = args.xc_functional
    force_eval["DFT-SCF-EPS_SCF"] = args.eps_scf
    force_eval["DFT-SCF-ADDED_MOS"] = args.added_mos
    force_eval["DFT-SCF-SMEAR"] = args.smear
    force_eval["DFT-SCF-SMEAR-METHOD"] = args.smear_method
    force_eval["DFT-SCF-SMEAR-ELECTRONIC_TEMPERATURE"] = args.electronic_temp
    force_eval["DFT-SCF-SMEAR-WINDOW_SIZE"] = args.window_size
    force_eval["DFT-SCF-DIAGONALIZATION"] = args.diag
    force_eval["DFT-SCF-OT"] = args.ot

    motion["MD-STEPS"] = args.md_steps
    motion["MD-ENSEMBLE"] = args.ensemble
    motion["MD-TEMPERATURE"] = args.temperature
    motion["PRINT-TRAJECTORY-FORMAT"] = args.traj_format

    task = md_run(args.file)
    task.md(directory=args.directory, runopt="genrun", force_eval=force_eval, motion=motion)
    task.analysis(directory=args.directory)

