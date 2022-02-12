#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.static import static_run

"""
usage:
"""

inputmopdos = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-static",
            help="directory for the static running")

    parser.add_argument("-f", "--file", type=str,
            help="the xyz file")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    parser.add_argument("--fileout", type=str, default="molecularpdos",
            help="prefix of output files containging molecular PDOS")

    parser.add_argument("--ngauss", type=int, default=0,
            help="gaussian broadening type")

    parser.add_argument("--degauss", type=float, default=0.001,
            help="gaussian broadening in Ry")

    parser.add_argument("--emin", type=str, default='default',
            help="min energy for PDOS")

    parser.add_argument("--emax", type=str, default='default',
            help="max energy for PDOS")

    parser.add_argument("--deltae", type=str, default='default',
            help="DeltaE: energy grid step (eV)")

    # -----------------------------------------------------------------
    #                       run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="molecularpdos",
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
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.molecularpdos(directory=args.directory, runopt=args.runopt, auto=args.auto)
