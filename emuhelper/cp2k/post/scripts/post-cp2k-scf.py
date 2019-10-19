#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from emuhelper.cp2k.post.scf import scf_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of static SCF running", type=str, default="tmp-cp2k-static")
    parser.add_argument("-f", "--file", help="output of scf running", type=str, default="static-scf.out")

    args = parser.parse_args()
    
    os.chdir(args.directory)
    task = scf_post(output=args.file)
    task.export()
    os.chdir("../")
