#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

from emuhelper.qe.static import static_run

"""
usage qe-converge-ecutwfc.py xxx.py emin emax step
"""

if __name__ == "__main__":
    task = static_run(sys.argv[1])
    task.converge_ecutwfc(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
