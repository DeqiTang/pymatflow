#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import sys

from emuhelper.qe.opt import opt_run

"""
usage: qe-opt.py xxx.xyz
"""

if __name__ == "__main__":
    task = opt_run(sys.argv[1])
    task.gen_input()
    task.run()
