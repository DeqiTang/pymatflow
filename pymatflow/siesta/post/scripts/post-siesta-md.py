#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.siesta.post.md import md_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of molecular dynamics running", type=str, default="tmp-siesta-md")
    parser.add_argument("-f", "--file", help="output of molecular dynamics running", type=str, default="molecular-dynamics.out")

    args = parser.parse_args()
    
    os.chdir(args.directory)
    task = md_post(outputfile=args.file)
    task.export()
    os.chdir("../")
