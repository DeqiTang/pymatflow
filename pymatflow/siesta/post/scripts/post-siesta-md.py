#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.siesta.post.md import MdPost


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of molecular dynamics running", type=str, default="tmp-siesta-md")
    parser.add_argument("-f", "--file", help="output of molecular dynamics running", type=str, default="molecular-dynamics.out")
    parser.add_argument("--view-traj", type=str, default="yes",
            choices=["yes", "no"],
            help="whether view the trajectory.")

    args = parser.parse_args()
    
    os.chdir(args.directory)
    task = MdPost(outputfile=args.file)
    task.export()
    if args.view_traj == "yes":
        task.view_trajectory()
    os.chdir("../")
