#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

from pymatflow.abinit.phonopy import phonopy

"""
usage:
    abinit-phonopy.py xxx.xyz
"""

if __name__ == "__main__":
    task = phonopy(sys.argv[1])
    task.gen_input()
    task.run()
    task.analysis()
