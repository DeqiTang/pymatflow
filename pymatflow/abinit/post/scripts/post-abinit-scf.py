#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.abinit.post.scf import scf_post



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously static running directory", type=str, default="tmp-abinit-static")
    parser.add_argument("--scfout", help="output file of static calculation", type=str, default="static-scf.out")

    args = parser.parse_args()

    os.chdir(args.directory)
    task = scf_post(args.scfout)
    task.export()
    os.chdir("../")
