#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.static import static_run

"""
usage:
    qe-scf.py -f xxx.xyz -k '2 2 2 0 0 0' --ecutwfc 100
"""


control_params = {}
system_params = {}
electrons_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("-d", "--directory", help="directory of the calculation", type=str, default="tmp-qe-static")
    parser.add_argument("--ecutwfc", help="ecutwfc, default value: 100 Ry", type=int, default=100)
    parser.add_argument("-k", "--kpoints", help="set kpoints like '3 3 3 0 0 0'", type=str, default="3 3 3 0 0 0")
 
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()
    xyzfile = args.file
    directory = args.directory
    system_params["ecutwfc"] = args.ecutwfc
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]
    

    task = static_run(xyzfile)
    task.scf(directory=directory, runopt="genrun", control=control_params, system=system_params, electrons=electrons_params, kpoints_mp=kpoints_mp)
