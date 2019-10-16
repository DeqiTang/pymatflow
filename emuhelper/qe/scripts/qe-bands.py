#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.static import static_run

"""
usage:
    qe-bands.py -f xxx.xyz
"""

control_params = {}
system_params = {}
electrons_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("-k", "--kptopt", help="kpoints schem option", type=str, default="automatic")
    parser.add_argument("--kpointsmp", help="the automatic schem kpoints", type=str, default="4 4 4 0 0 0")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file
    kptopt = args.kptopt
    kpoints_mp = [int(args.kpointsmp.split()[i]) for i in range(6)]


    task = static_run(xyzfile)
    task.bands(kptopt=kptopt, control=control_params, system=system_params, electrons=electrons_params, kpoints_mp=kpoints_mp)
