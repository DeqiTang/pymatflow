#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

from emuhelper.siesta.md import md_run

"""
usage:
    siesta-md.py xxx.xyz
"""


if __name__ == "__main__":
    task = md_run(sys.argv[1])
    task.gen_input()
    task.run()
    task.analysis()
