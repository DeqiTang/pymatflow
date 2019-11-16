#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.siesta.post.converge import converge_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of converge test running", type=str, default="tmp-siesta-cutoff")
    parser.add_argument("-c", "--converge", help="type of converge test(MeshCutoff)", type=str, default="MeshCutoff")
    args = parser.parse_args()
    
    task = converge_post()
    task.postprocess(directory=args.directory, converge=args.converge)
