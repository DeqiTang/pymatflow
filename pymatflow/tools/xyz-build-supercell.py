#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from pymatflow.base.xyz import base_xyz

"""
usage:
    build_supercell input.xyz output.xyz n1 n2 n3
    n1, n2, n3 are the three repeated number in
    the three direct of three basis vectors
"""


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, 
            help="input xyz file")

    parser.add_argument("-o", "--output", type=str,
            help="output xyz file")

    parser.add_argument("-n", "--supern", nargs="+", type=int,
            help="bulid supern:[int, int, int] supercell")

    args = parser.parse_args()

    xyz = base_xyz()
    xyz.get_xyz(args.input)
    xyz.build_supercell(args.supern)
    xyz.to_xyz(args.output)
