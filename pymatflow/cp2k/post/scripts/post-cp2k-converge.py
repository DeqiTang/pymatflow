#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.cp2k.post.converge import converge_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of converge test running", type=str, default="tmp-cp2k-cutoff")
    parser.add_argument("-c", "--converge", help="type of converge test(cutoff, rel_cutoff)", type=str, default="cutoff")
    args = parser.parse_args()
    
    task = converge_post()
    task.postprocess(directory=args.directory, converge=args.converge)
