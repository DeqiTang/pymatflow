#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse 
from emuhelper.qe.static import static_run

"""
usage:
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of the static running", type=str, default="tmp-qe-static")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--runopt", help="gen, run, or genrun", type=str, default="genrun")
    parser.add_argument("--plot-num", type=int, default=0,
            choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17, 18, 19, 20, 21],
            help="type of analysis stored in the filplot file for later plot")
    parser.add_argument("--iflag", type=int,
            default=3,
            choices=[0, 1, 2, 3, 4],
            help="dimension of the plot")
    parser.add_argument("--output-format", type=int, default=5,
            help="output file format for visualization")
 
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()
    xyzfile = args.file

 
    task = static_run(xyzfile)
    task.pp(directory=args.directory, runopt=args.runopt, plot_num=args.plot_num, iflag=args.iflag, output_format=args.output_format)
