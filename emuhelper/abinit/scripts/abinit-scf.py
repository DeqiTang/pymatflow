#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import argparse
from emuhelper.abinit.static import static_run

"""
usage:
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", type=str, default="tmp-abinit-static",
            help="Directory for the static running.")
    parser.add_argument("-f", "--file", type=str,
            help="The xyz file name.")
    parser.add_argument("--runopt", type=str, default="genrun",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")
    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")
    parser.add_argument("--ecut", type=int, default=15,
            help="Kinetic energy cutoff for wave functions in unit of Hartree, default value: 15 Hartree")
    
    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
    # ==========================================================   
    args = parser.parse_args()

    task = static_run(args.file)
    task.scf(directory=args.directory, mpi=args.mpi, runopt=args.runopt)
