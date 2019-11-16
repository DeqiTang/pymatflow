#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.qe.post.scf import scf_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of static scf running", type=str, default="tmp-qe-static")
    parser.add_argument("-f", "--file", help="output of static scf running", type=str, default="static-scf.out")

    args = parser.parse_args()
    
    os.chdir(args.directory)
    task = scf_post(output=args.file)
    task.export()
    os.chdir("../")
