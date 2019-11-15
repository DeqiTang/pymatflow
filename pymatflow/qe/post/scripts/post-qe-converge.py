#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from emuhelper.qe.post.converge import converge_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of static scf running", type=str, default="tmp-qe-ecutwfc")
    parser.add_argument("-c", "--converge", help="type of converge test(ecutwfc, ecutrho, kpoints)", type=str)
    args = parser.parse_args()
    
    task = converge_post()
    task.postprocess(directory=args.directory, converge=args.converge)
