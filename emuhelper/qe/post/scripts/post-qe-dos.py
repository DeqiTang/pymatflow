#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

from emuhelper.qe.post.dos import qe_dos

"""
usage: post-qe-dos.py xxx.dos
"""

if __name__ == "__main__":
    task = qe_dos(sys.argv[1])
    task.plot_dos()
