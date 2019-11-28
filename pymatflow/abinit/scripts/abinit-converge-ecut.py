#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import argparse 

from pymatflow.abinit.static import static_run

"""
usage: abinit-converge-ecut.py xxx.xyz emin emax step
"""

electrons_params = {}
kpoints_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-abinit-ecut",
            help="Directory for the ecut converge running.")

    parser.add_argument("-f", "--file", type=str,
            help="The xyz file name.")

    parser.add_argument("--range", nargs="+", type=int,
            help="ecut test range")

    parser.add_argument("--runopt", type=str, default="genrun",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--ixc", type=int, default=11,
            choices=[1, 2, 3 ,4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 40, 41, 42],
            help="type of exchage-correlation functional. for more information, refer to https://docs.abinit.org/variables/basic/#ixc")

    parser.add_argument("--kptopt", type=int, default=1,
            help="Kpoints Generation scheme option: 0, 1, 2, 3, 4 or a negative value. for more information, refer to https://docs.abinit.org/variables/basic/#kptopt")

    parser.add_argument("--ngkpt", nargs="+", type=int,
            default=[1, 1, 1],
            help="number of grid points for kpoints generation. for more information, refer to https://docs.abinit.org/variables/basic/#ngkpt")
    
    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
    # ==========================================================   
    args = parser.parse_args()

    electrons_params["ixc"] = args.ixc

    kpoints_params["kptopt"] = args.kptopt
    kpoints_params["ngkpt"] = args.ngkpt

    
    task = static_run(args.file)
    task.converge_ecut(args.range[0], args.range[1], args.range[2], directory=args.directory, mpi=args.mpi, runopt=args.runopt, electrons=electrons_params, kpoints=kpoints_params)
