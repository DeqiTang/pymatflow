#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import sys

from emuhelper.cp2k.static import static_run

"""
usage: cp2k-bands.py xxx.xyz
"""

if __name__ == "__main__":
    task = static_run(sys.argv[1])
    task.print_bands()
    task.gen_input()
    task.run()
