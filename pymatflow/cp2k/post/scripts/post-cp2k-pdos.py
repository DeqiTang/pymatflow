#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.cp2k.post.pdos import pdos_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of static SCF running", type=str, default="tmp-cp2k-static")
    parser.add_argument("-f", "--files", help="output of pdos", nargs="+", type=str)
    parser.add_argument("--sigma", help="Sigma", type=float, default=0.01)
    parser.add_argument("--step", help="new energy step", type=float, default=0.001)

    args = parser.parse_args()
    
    os.chdir(args.directory)
    task = pdos_post(pdos_files=args.files, step=args.step, sigma=args.sigma)
    task.export()
    os.chdir("../")
