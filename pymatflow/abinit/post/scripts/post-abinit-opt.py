#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.abinit.post.opt import opt_post



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously optimization running directory", type=str, default="tmp-abinit-opt")
    parser.add_argument("--optout", help="output file of optimization calculation", type=str, default="geometric-optimization.out")

    args = parser.parse_args()

    os.chdir(args.directory)
    task = opt_post(args.optout)
    task.export()
    os.chdir("../")
