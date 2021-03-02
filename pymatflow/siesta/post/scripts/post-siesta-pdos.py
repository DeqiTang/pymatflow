#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

from pymatflow.siesta.post.pdos import Pdos

"""
usage: post-siesta-pdos.py xxx.PDOS.xml
"""

if __name__ == "__main__":
    task = Pdos(sys.argv[1])
    task.projected_to_element()
