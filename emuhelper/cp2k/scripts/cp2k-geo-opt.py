#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import os
import shutil

import matplotlib.pyplot as plt

from emuhelper.cp2k.opt import opt_run

"""
Usage:
    cp2k-geo-opt.py xxx.xyz
    xxx.xyz is the input structure file

"""

if __name__ == "__main__":
    task = opt_run(sys.argv[1])
    task.geo_opt(directory="tmp-cp2k-geo-opt", runopt="genrun")
    task.analysis(directory="tmp-cp2k-geo-opt")



