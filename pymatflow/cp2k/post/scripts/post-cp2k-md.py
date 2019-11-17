#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.cp2k.post.md import md_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of Molecular Dynamics running", type=str, default="tmp-cp2k-md")
    parser.add_argument("-f", "--file", help="output of opt running", type=str, default="md.out")

    args = parser.parse_args()
    
    os.chdir(args.directory)
    task = md_post(output=args.file, run_type='MD')
    task.export()
    os.chdir("../")
