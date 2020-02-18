#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.static import static_run
from pymatflow.remote.server import server_handle

"""
usage:
"""

inputmopdos = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory for the static running", type=str, default="tmp-qe-static")
    parser.add_argument("--mpi", help="MPI commadn", type=str, default="")
    parser.add_argument("-f", "--file", help="the xyz file", type=str)
    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")
    parser.add_argument("--fileout", help="prefix of output files containging molecular PDOS", type=str, default="molecularpdos")
    parser.add_argument("--ngauss", help="gaussian broadening type", type=int, default=0)
    parser.add_argument("--degauss", help="gaussian broadening in Ry", type=float, default=0.001)
    parser.add_argument("--emin", help="min energy for PDOS", type=str, default='default')
    parser.add_argument("--emax", help="max energy for PDOS", type=str, default='default')
    parser.add_argument("--deltae", help="DeltaE: energy grid step (eV)", type=str, default='default')

    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")
    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")
    parser.add_argument("--jobname", type=str, default="pwscf-scf",
            help="jobname on the pbs server")
    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")
    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")



    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file

    if args.emin == 'default':
        emin = args.emin
    else:
        emin = float(args.emin)
    if args.emax == 'default':
        emax = args.emax
    else:
        emax = float(args.emax)
    if args.deltae == 'default':
        deltae = args.deltae
    else:
        deltae = float(args.deltae)

    inputmopdos["fileout"] = args.fileout
    inputmopdos["ngauss"] = args.ngauss
    inputmopdos["degauss"] = args.degauss
    inputmopdos["emin"] = emin
    inputmopdos["emax"] = emax
    inputmopdos["deltae"] = deltae

    task = static_run()
    task.get_xyz(args.file)

    task.set_molecularpdos(inputmopdos=inputmopdos)
    task.molecularpdos(directory=args.directory, runopt=args.runopt, mpi=args.mpi, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)

    server_handle(auto=args.auto, directory=args.directory, jobfilebase="static-molecularpdos", server=args.server)
