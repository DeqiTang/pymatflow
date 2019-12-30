#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from pymatflow.qe.static import static_run

"""
usage:
    qe-turbo-davidson-spectrum.py -f xxx.xyz -k '2 2 2 0 0 0' --ecutwfc 100
"""



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="the xyz file", type=str)
    parser.add_argument("-d", "--directory", help="directory of the calculation", type=str, default="tmp-qe-static")
    parser.add_argument("--runopt", help="gen, run, or genrun", type=str, default="genrun")
 
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()
    directory = args.directory
    xyzfile = args.file
    
    task = static_run()
    task.get_xyz(args.file)
    task.turbo_davidson(directory=directory, runopt=args.runopt)
