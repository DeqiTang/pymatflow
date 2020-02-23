#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.static import static_run


"""
usage:
    qe-pdos.py -f xxx.xyz
"""

projwfc_input = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-static",
            help="directory for the static running")

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command")

    parser.add_argument("-f", "--file", type=str,
            help="the xyz file")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    parser.add_argument("--filpdos", type=str, default="projwfc",
            help="output projected dos file name")

    parser.add_argument("--ngauss", type=str, default="default",
            help="gaussian broadening type")

    parser.add_argument("--degauss", type=str, default='default',
            help="gaussian broadening")

    parser.add_argument("--emin", type=str, default='default',
            help="min energy for DOS")

    parser.add_argument("--emax", type=str, default='default',
            help="max energy for DOS")

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

    parser.add_argument("--jobname", type=str, default="projwfc-pdos",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

    if args.ngauss == 'default':
        ngauss = args.ngauss
    else:
        ngauss = int(args.ngauss)
    if args.degauss == 'default':
        degauss = args.degauss
    else:
        degauss = float(args.degauss)
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

    projwfc_input["filpdos"] = args.filpdos
    projwfc_input["ngauss"] = ngauss
    projwfc_input["degauss"] = degauss
    projwfc_input["emin"] = emin
    projwfc_input["emax"] = emax
    projwfc_input["deltae"] = deltae

    task = static_run()
    task.get_xyz(args.file)
    task.set_projwfc(projwfc_input=projwfc_input)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.projwfc(directory=args.directory, runopt=args.runopt, auto=args.auto)
