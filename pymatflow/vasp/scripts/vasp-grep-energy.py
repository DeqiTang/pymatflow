#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse


"""
usage:
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, default="./OUTCAR",
            help="the OUTCAR file")

    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
    # ==========================================================
    args = parser.parse_args()

    os.system("cat %s | grep \"entropy=\"" % (args.input))