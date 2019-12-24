#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from pymatflow.qe.post.pdos import pdos_post

"""
usage: post-qe-pdos.py 
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously projwfc running directory", type=str, default="tmp-qe-static")
    parser.add_argument("--filpdos", help="filpdos name used in projwfc running", type=str, default="projwfc")

    parser.add_argument("--plotrange", type=float, nargs="+",
            default=[0, 1.0],
            help="plot range (in percentage), like --plotrange 0.1 0.9")

    parser.add_argument("--atomtoproj", type=int, nargs="+",
            default=[],
            help="atom to projection in atom projected dos. atom number starts with 1.")

    parser.add_argument("--fontsize", type=int, default=10,
            help="fontsize for the plot.")

    args = parser.parse_args()

    directory = args.directory
    filpdos = args.filpdos

    task = pdos_post()
    task.get_data(directory=directory, filpdos=filpdos)
    task.export(directory=directory, plotrange=args.plotrange, atomtoproj=args.atomtoproj, fontsize=args.fontsize)
