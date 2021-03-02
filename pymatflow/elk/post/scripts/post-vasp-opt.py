#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.vasp.post.opt import opt_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of relax running", type=str, default="tmp-vasp-optimization")
    parser.add_argument("-f", "--file", help="output of opt running", type=str, default="OUTCAR")
    parser.add_argument("--view-traj", type=str, default="yes",
            choices=["yes", "no"],
            help="whether to view the trajectory.")

    args = parser.parse_args()

    os.chdir(args.directory)
    task = opt_post(output=args.file)
    task.export()
    if args.view_traj == "yes":
        task.view_trajectory()
    os.chdir("../")
