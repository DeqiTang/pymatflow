#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

from emuhelper.base.xyz import base_xyz

"""
usage:
    build_supercell input.xyz output.xyz n1 n2 n3
    n1, n2, n3 are the three repeated number in
    the three direct of three basis vectors
"""


if __name__ == '__main__':
    xyz = base_xyz(sys.argv[1])
    supern = [int(sys.argv[i]) for i in range(3, 6)]
    xyz.build_supercell(supern)
    xyz.to_xyz(sys.argv[2])
