#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import os
import shutil
import getopt

import matplotlib.pyplot as plt

from emuhelper.cp2k.opt import opt_run

"""
Usage:
    cp2k-opt.py xxx.xyz -h -t xxx ( or --type xxx )
    xxx.xyz is the input structure file

    if you specify no options and args, the script
    will conduct GEO_OPT by default.
"""

def parse_args(argv, task):
    """
    task: an object of opt_run
    """
    try:
        opts, args = getopt.getopt(argv, "t:h", ["type="])
    except getopt.GetoptError:
        print("cp2k-opt.py -h\n")
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print("cp2k-opt.py input.xyz --options\n")
            print("if you specify no options and args\n")
            print("the scripts will conduct GEO_OPT\n")
            print("by default\n")
            sys.exit(0)
        elif opt == '-t' or opt == '--type':
            task.set_run_type(arg)

if __name__ == "__main__":
    task = opt_run(sys.argv[1])
    #task.set_run_type("GEO_OPT")
    #task.set_run_type("CELL_OPT")
    parse_args(sys.argv[2:], task)
    task.gen_input()
    task.run()
    task.analysis()



