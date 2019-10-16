#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.static import static_run

"""
usage:
    qe-pdos.py -f xxx.xyz
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--filpdos", help="output projected dos file name", type=str, default="output.pdos")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file
    filpdos = args.filpdos

    task = static_run(xyzfile)
    task.projwfc(filpdos=filpdos)
