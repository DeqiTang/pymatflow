#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.cp2k.post.opt import opt_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of CELL_OPT running", type=str, default="tmp-cp2k-cell-opt")
    parser.add_argument("-f", "--file", help="output of opt running", type=str, default="cell-opt.out")
    #parser.add_argument("-o", "--output", help="converted to xyz file", type=str, default="relaxed.xyz")

    args = parser.parse_args()
    
    os.chdir(args.directory)
    task = opt_post(output=args.file, run_type='CELL_OPT')
    task.export()
    os.chdir("../")
