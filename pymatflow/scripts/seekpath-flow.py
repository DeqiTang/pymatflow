#!/usr/bin/env python
#import os
import sys
#import numpy as np
#from pymatflow.base.xyz import BaseXyz
from pymatflow.cmd.structflow import read_structure
from pymatflow.cmd.structflow import write_structure
from pymatflow.base.element import element
#import seekpath
#import spglib
import argparse

"""
Reference:
    https://atztogo.github.io/spglib/python-spglib.html
    https://seekpath.readthedocs.io/en/latest/index.html
Warning:
    the result is not guaranteed to be correct
"""



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="driver", title="subcommands", description="choose one and only one sub command")

    # ------------------------
    # info
    # ------------------------
    subparser = subparsers.add_parser("info", help="print out the info of the structure")

    gp = subparser.add_argument_group(title="Input/Output")

    gp.add_argument("-i", "--input", type=str,
            help="the input structure file")

    gp = subparser.add_argument_group(title="seekpath parameters")

    gp.add_argument("--with-time-reversal", type=str, default="true",
            choices=["True", "true", "False", "false"],
            help="if time-reversal symmetry is present or not, default is true")

    gp.add_argument("--symprec", type=float, default=1.0e-9,
            help="symmetry precision")

    # --------------------------
    # std
    # --------------------------
    subparser = subparsers.add_parser("std", help="standardize the structure using seekpath")

    gp = subparser.add_argument_group(title="Input/Output")

    gp.add_argument("-i", "--input", type=str,
        help="the input structure file")

    gp.add_argument("-o", "--output", type=str,
        help="the output standarized structure file")

    gp = subparser.add_argument_group(title="seekpath parameters")

    gp.add_argument("--with-time-reversal", type=str, default="true",
            choices=["True", "true", "False", "false"],
            help="if time-reversal symmetry is present or not, default is true")

    gp.add_argument("--symprec", type=float, default=1.0e-9,
            help="symmetry precision")
        
    # --------------------------
    # kpath
    # --------------------------
    subparser = subparsers.add_parser("kpath", help="get the kpath from seekpath")

    gp = subparser.add_argument_group(title="Input/Output")

    gp.add_argument("-i", "--input", type=str,
        help="the input structure file")

    gp.add_argument("-o", "--output", type=str, default="kpath-from-seekpath.txt",
            help="the output kpoitns file")

    gp = subparser.add_argument_group(title="seekpath parameters")

    gp.add_argument("--with-time-reversal", type=str, default="true",
            choices=["True", "true", "False", "false"],
            help="if time-reversal symmetry is present or not, default is true")

    gp.add_argument("--symprec", type=float, default=1.0e-9,
            help="symmetry precision")    

    gp.add_argument("--join", type=int, default=15,
            help="default number of kpoint to connect the connected high symmetry k point")

    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()

    # if no argument passed to matflow
    if len(sys.argv) == 1:
        # display help message when no args provided
        parser.print_help()
        sys.exit(1)

    if args.driver == "info":
        from pymatflow.third.seekpath import get_spacegroup_and_kpath
        a = read_structure(filepath=args.input)
        print("===========================================\n")
        print("calculated using spglib\n")
        print("===========================================\n")
        print("spacegroup is : %s\n" % get_spacegroup_and_kpath(structure=a, symprec=args.symprec)[0])
        print("Warning:\n")
        print("the result is not guaranteed to be correct\n")
        print("-------------------------------------------\n")
        print("suggested k path calculated using seekpath:\n")
        with_time_reversal = True if args.with_time_reversal in ["True", "true"] else False
        print(get_spacegroup_and_kpath(structure=a, with_time_reversal=with_time_reversal, symprec=args.symprec)[1])
    elif args.driver == "std":
        from pymatflow.third.seekpath import seekpath_std_structure
        with_time_reversal = True if args.with_time_reversal in ["True", "true"] else False
        std_structure = seekpath_std_structure(
            structure=read_structure(filepath=args.input),
            with_time_reversal=with_time_reversal,
            symprec=args.symprec
        )
        write_structure(std_structure, args.output)
    elif args.driver == "kpath":
        from pymatflow.third.seekpath import seekpath_output_kpath
        structure = read_structure(filepath=args.input)
        with_time_reversal = True if args.with_time_reversal in ["True", "true"] else False
        seekpath_output_kpath(
            structure=structure,
            output=args.output,
            with_time_reversal=with_time_reversal,
            symprec=args.symprec,
            join=args.join
        )
    else:
        pass