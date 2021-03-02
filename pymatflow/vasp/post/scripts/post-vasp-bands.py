#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

import numpy as np
import matplotlib.pyplot as plt

from pymatflow.vasp.post.bands import post_bands


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-vasp-static",
            help="directory of the static calculation")

    parser.add_argument("--kpath", type=str, nargs="+", default=None,
            help="manual input kpath for band structure calculation")

    parser.add_argument("--kpath-file", type=str, default="kpath-from-seekpath.txt",
            help="file to read the kpath for band structure calculation")

    parser.add_argument("--bands-p4vasp", type=str, default=None,
            help="bands data exported by p4vasp in gnuplot format")

    parser.add_argument("--bandrange", type=float, nargs="+",
            default=[0, 1.0],
            help="band range to plot. in percentage")

    parser.add_argument("--engine", type=str, default="matplotlib",
            choices=["matplotlib", "gnuplot"],
            help="plot engine, matplotlib or gnuplot")

    parser.add_argument("--output-dir", type=str, default="./",
           help="directory to put the analysis result files")


    args = parser.parse_args()

    bands = post_bands()
    bands.get_vasprun(os.path.join(args.directory, "vasprun.xml"))
    bands.get_kpath(kpath_manual=args.kpath, kpath_file=args.kpath_file)
    bands.export(directory=args.directory, option=args.engine, bandrange=args.bandrange)
