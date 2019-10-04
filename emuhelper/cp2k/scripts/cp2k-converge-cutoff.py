#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

from emuhelper.cp2k.static import static_run

"""
usage: cp2k-converge-cutoff.py xxx.xyz emin emax step rel_cutoff
"""

if __name__ == "__main__":
    task = static_run(sys.argv[1])
    task.converge_cutoff(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))

