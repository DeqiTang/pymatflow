#!/usr/bin/env python

import argparse

from pymatflow.base.xyz import BaseXyz
from pymatflow.cmd.structflow import read_structure

from pymatflow.third.seekpath import seekpath_get_kpath

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    parser.add_argument("-o", "--output", type=str, default="kpath-from-seekpath.txt",
            help="the output kpoitns file")

    parser.add_argument("--with-time-reversal", type=str, default="true",
            choices=["True", "true", "False", "false"],
            help="if time-reversal symmetry is present or not, default is true")

    parser.add_argument("--symprec", type=float, default=1.0e-9,
            help="symmetry precision")

    parser.add_argument("--join", type=int, default=15,
            help="default number of kpoint to connect the connected high symmetry k point")

    # ===============================================================================
    args = parser.parse_args()
    
    structure = read_structure(filepath=args.input)
    
    with_time_reversal = True if args.with_time_reversal in ["True", "true"] else False
    seekpath_get_kpath(
        structure=structure,
        output=args.output,
        with_time_reversal=with_time_reversal,
        symprec=args.symprec,
        join=args.join
    )


if __name__ == "__main__":
    main()
