#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from emuhelper.qe.post.opt import opt_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of relax running", type=str, default="tmp-qe-relax")
    parser.add_argument("-f", "--file", help="output of opt running", type=str, default="relax.out")
    parser.add_argument("-o", "--output", help="converted to xyz file", type=str, default="relaxed.xyz")

    args = parser.parse_args()
    
    os.chdir(args.directory)
    task = opt_post(output=args.file, run_type='relax')
    task.export()
    os.chdir("../")
