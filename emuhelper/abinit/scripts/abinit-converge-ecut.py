#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

from emuhelper.abinit.static import static_run

"""
usage: abinit-converge-ecut.py xxx.xyz emin emax step
"""

if __name__ == "__main__":
    task = static_run(sys.argv[1])
    task.converge_ecut(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
