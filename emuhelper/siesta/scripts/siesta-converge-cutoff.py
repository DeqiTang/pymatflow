#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

from emuhelper.siesta.static import static_run

"""
usage:
    siesta-converge-cutoff.py xxx.xyz emin emax step
"""

if __name__ == "__main__":
    task = static_run(sys.argv[1])
    task.converge_cutoff(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))

