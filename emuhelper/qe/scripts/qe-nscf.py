#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
from emuhelper.qe.static import static_run

"""
usage:
    qe-single-point.py xxx.xyz
"""

if __name__ == "__main__":
    task = static_run(sys.argv[1])
    task.nscf()
    task.run_nscf()
