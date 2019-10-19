#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.cp2k.static import static_run

"""
usage: cp2k-converge-rel-cutoff.py xxx.xyz emin emax step cutoff
"""

force_eval = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of the calculation", type=str, default="tmp-cp2k-static")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--xc-functional", help="XC_FUNCTIONAL type", type=str, default="PBE")
    parser.add_argument("--range", help="REL_CUTOFF range emin emax step", nargs='+', type=int)
    parser.add_argument("--cutoff", help="REL_CUTOFF, default value: 60 Ry", type=int , default=60)
    parser.add_argument("-k", "--kpoints", help="set kpoints like '3 3 3 0 0 0'", type=str, default="3 3 3 0 0 0")
    parser.add_argument("--added-mos", help="ADDED_MOS for SCF", type=int, default=0)
    parser.add_argument("--smear", help="smearing type", type=str, default="FERMI_DIRAC")
    parser.add_argument("--electronic-temp", help="ELECTRON_TEMPERATURE for FERMI_DIRAC SMEAR", type=float, default=300)
    parser.add_argument("--window-size", help="Size of the energy window centred at the Fermi level for ENERGY_WINDOW type smearing", type=float, default=0)
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()
    directory = args.directory
    xyzfile = args.file
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]
    force_eval["XC_FUNCTIONAL"] = args.xc_functional
    force_eval["ADDED_MOS"] = args.added_mos
    force_eval["SMEAR"] = args.smear
    force_eval["ELECTRONIC_TEMPERATURE"] = args.electronic_temp
    force_eval["WINDOW_SIZE"] = args.window_size

    task = static_run(xyzfile)
    task.converge_rel_cutoff(args.range[0], args.range[1], args.range[2], cutoff=args.cutoff, force_eval=force_eval, runopt="genrun")
