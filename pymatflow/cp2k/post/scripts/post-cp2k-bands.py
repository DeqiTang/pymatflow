#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.cp2k.post.bands import bands_post

"""
usage:
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously static running directory", type=str, default="tmp-cp2k-static")
    parser.add_argument("--bandsfile", help="output band file", type=str, default="bands.bs")
    parser.add_argument("--option", help="choosing gnuplot or matplotlib to do the band plot", type=str, default="gnuplot")

    args = parser.parse_args()


    os.chdir(args.directory)
    task = bands_post()
    task.plot_band(bandsfile=args.bandsfile, option=args.option)
