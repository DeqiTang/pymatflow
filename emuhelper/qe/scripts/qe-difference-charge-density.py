#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
from emuhelper.qe.static import static_run

"""
usage:
    qe-difference-charge-density.py xxx.xyz
"""

if __name__ == "__main__":
    task = static_run(sys.argv[1])
    task.difference_charge_density()
