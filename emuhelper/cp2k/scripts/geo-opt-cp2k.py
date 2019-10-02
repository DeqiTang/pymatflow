#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import os
import shutil

import matplotlib.pyplot as plt

from emuhelper.cp2k.opt import opt_run

"""
Usage:
    python geo_opt_cp2k.py xxx.xyz 
    xxx.xyz is the input structure file

    make sure the xyz structure file 
    is in the directory.
"""


if __name__ == "__main__":
    opt = opt_run(sys.argv[1])
    opt.gen_input()
    opt.run()
    opt.analysis()



