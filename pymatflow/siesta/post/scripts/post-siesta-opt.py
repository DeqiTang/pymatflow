#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.siesta.post.opt import opt_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of geometric optimization running", type=str, default="tmp-siesta-opt")
    parser.add_argument("-f", "--file", help="output of opt running", type=str, default="geometric-optimization.out")

    args = parser.parse_args()
    
    os.chdir(args.directory)
    task = opt_post(outputfile=args.file)
    task.export()
    os.chdir("../")
