#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.qe.post.converge import ConvergePost


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of static scf running", type=str, default="tmp-qe-ecutwfc")
    parser.add_argument("-c", "--converge", help="type of converge test(ecutwfc, ecutrho, kpoints)", type=str)
    args = parser.parse_args()
    
    task = ConvergePost()
    task.postprocess(directory=args.directory, converge=args.converge)
