#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

from emuhelper.abinit.opt import opt_run

"""
usage: abinit-opt.py xxx.xyz
"""

if __name__ == "__main__":
    task = opt_run(sys.argv[1])
    task.gen_input()
    task.run()

