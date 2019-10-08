#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

from emuhelper.abinit.static import static_run

"""
usage:
    abinit-dft-u.py xxx.xyz
"""

if __name__ == '__main__':
    task = static_run(sys.argv[1])
    task.dft_plus_u()
    task.gen_input()
    task.run()
