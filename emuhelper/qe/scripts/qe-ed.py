#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse 
from emuhelper.qe.static import static_run

"""
usage:
    qe-difference-charge-density.py xxx.xyz
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("-d", "--directory", help="directory of the calculation", type=str, default="tmp-qe-static")
 
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()
    xyzfile = args.file
    directory = args.directory
 
    task = static_run(xyzfile)
    task.electron_density()
