#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.post.bands import bands_post

"""
usage:
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously projwfc running directory", type=str, default="tmp-qe-static")
    parser.add_argument("--pwxbandsin", help="input file for the pw.x band calculation", type=str, default="static-bands.in")
    parser.add_argument("--bandsxout", help="outputfile of the band,x calculation", type=str, default="bands.out")
    parser.add_argument("--option", help="choosing gnuplot or matplotlib to do the band plot", type=str, default="gnuplot")

    parser.add_argument("--bandrange", type=float, nargs="+",
            default=[0, 1.0],
            help="band range to plot. in percentage")

    args = parser.parse_args()


    os.chdir(args.directory)
    task = bands_post(pwxbandsin=args.pwxbandsin, bandsxout=args.bandsxout)
    task.plot_band(option=args.option, bandrange=args.bandrange)
