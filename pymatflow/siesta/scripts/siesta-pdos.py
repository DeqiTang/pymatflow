#!/usr/bin/env python
# _*_ coding: utf-8

import sys

from pymatflow.siesta.static import static_run

"""
usage: siesta-pdos.py xxx.xyz
"""

if __name__ == "__main__":
    task = static_run(sys.argv[1])
    task.properties.option = "pdos"
    #task.set_spin("polarized")
    task.gen_input()
    task.run()
    task.analysis()
