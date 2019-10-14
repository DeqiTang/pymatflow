#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.static import static_run

"""
usage qe-converge-ecutrho.py -f xxx.py --range emin emax step --ecutwfc xxx
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--range", help="ecutrho test range", nargs='+', type=int)
    parser.add_argument("--ecutwfc", help="ecutwfc value", type=int)
    parser.add_argument("-k", "--kpoints", help="set kpoints like '1 1 1 0 0 0'", type=str, default="1 1 1 0 0 0")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]
    
    task = static_run(xyzfile)
    task.converge_ecutrho(args.range[0], args.range[1], args.range[2], args.ecutwfc, runopt="genrun")
