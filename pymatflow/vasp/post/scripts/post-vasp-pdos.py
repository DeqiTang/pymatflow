#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

import numpy as np
import matplotlib.pyplot as plt

from pymatflow.vasp.post.pdos import post_pdos


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-vasp-static",
            help="directory of the static calculation")

    parser.add_argument("--plotrange", type=float, nargs="+",
            default=[0, 1.0],
            help="range to plot. in percentage")

    parser.add_argument("--engine", type=str, default="matplotlib",
            choices=["matplotlib", "gnuplot"],
            help="plot engine, matplotlib or gnuplot")


    args = parser.parse_args()

    pdos = post_pdos()
    pdos.get_vasprun(os.path.join(args.directory, "vasprun.xml"))
    pdos.export(directory=args.directory, option=args.engine, plotrange=args.plotrange)
