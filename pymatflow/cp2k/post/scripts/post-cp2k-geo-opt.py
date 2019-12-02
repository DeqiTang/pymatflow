#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.cp2k.post.opt import opt_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of GEO_OPT running", type=str, default="tmp-cp2k-geo-opt")
    parser.add_argument("-f", "--file", help="output of opt running", type=str, default="geo-opt.out")
    #parser.add_argument("-o", "--output", help="converted to xyz file", type=str, default="relaxed.xyz")
    parser.add_argument("--view-traj", type=str, default="yes",
            choices=["yes", "no"],
            help="whether view the trajectory.")

    args = parser.parse_args()
    
    os.chdir(args.directory)
    task = opt_post(output=args.file, run_type='GEO_OPT')
    task.export()
    if args.view_traj == "yes":
        task.view_trajectory()
    os.chdir("../")
