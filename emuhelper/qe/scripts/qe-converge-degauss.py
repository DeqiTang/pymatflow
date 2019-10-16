#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.static import static_run

"""
usage qe-converge-degauss.py -f xxx.xyz --range degauss_min degauss_max step
"""

control_params = {}
# do not set occupation, smearing and degauss through system_params
system_params = {}
electrons_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--range", help="degauss_min degauss_max step", nargs='+', type=float)
    parser.add_argument("-k", "--kpoints", help="set kpoints like '1 1 1 0 0 0'", type=str, default="1 1 1 0 0 0")
    parser.add_argument("--ecutwfc", help="better a previously converged ecutwfc", type=int, default=100)

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]
    system_params["ecutwfc"] = args.ecutwfc

    task = static_run(xyzfile)

    task.converge_degauss(round(args.range[0], 6), round(args.range[1], 6), round(args.range[2], 6), control=control_params, system=system_params, electrons=electrons_params, runopt="genrun")
