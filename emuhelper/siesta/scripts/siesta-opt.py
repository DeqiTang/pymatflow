#!/usr/bin/env python
 # _*_ coding: utf-8 _*_

import sys

from emuhelper.siesta.opt import opt_run

"""
usage:
   siesta-opt.py xxx.xyz
"""

if __name__ == "__main__":
    task = opt_run(sys.argv[1])
    task.gen_input()
    task.run()
    task.analysis()
