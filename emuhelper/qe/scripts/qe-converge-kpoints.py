#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

from emuhelper.qe.static import static_run

"""
usage qe-converge-kpoints.py xxx.xyz nk_min nk_max step
"""

if __name__ == "__main__":
    task = static_run(sys.argv[1])
    task.converge_kpoints(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
