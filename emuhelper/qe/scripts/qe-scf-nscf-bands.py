#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.static import static_run

"""
usage:
    qe-single-point.py xxx.xyz
"""

control_params = {}
system_params = {}
electrons_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--ecutwfc", help="ecutwfc, default value: 100 Ry", type=int, default=100)
    parser.add_argument("--scfkpoints", help="set scf kpoints like '3 3 3 0 0 0'", type=str, default="3 3 3 0 0 0")
    parser.add_argument("--nscfkpoints", help="set nscf kpoints like '3 3 3 0 0 0'", type=str, default="4 4 4 0 0 0")
    parser.add_argument("-k", "--kptopt", help="kpoints schem option", type=str, default="automatic")
    parser.add_argument("--bandskpoints", help="the automatic schem kpoints", type=str, default="6 6 6 0 0 0")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file
    system_params["ecutwfc"] = args.ecutwfc
    scf_kpoints_mp = [int(args.scfkpoints.split()[i]) for i in range(6)]
    nscf_kpoints_mp = [int(args.nscfkpoints.split()[i]) for i in range(6)]
    bands_kpoints_mp = [int(args.nscfkpoints.split()[i]) for i in range(6)]
    kptopt = args.kptopt

    task = static_run(xyzfile)
    task.scf(runopt="genrun", control=control_params, system=system_params, electrons=electrons_params, kpoints_mp=scf_kpoints_mp)    
    task.nscf(runopt='genrun', system=system_params, electrons=electrons_params, kpoints_mp=nscf_kpoints_mp)
    task.bands(kptopt=kptopt, control=control_params, system=system_params, electrons=electrons_params, kpoints_mp=bands_kpoints_mp)
