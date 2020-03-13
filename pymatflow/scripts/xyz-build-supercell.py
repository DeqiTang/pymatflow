#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from pymatflow.base.xyz import base_xyz
from pymatflow.structure.crystal import crystal

"""
usage:
    xyz-build-supercell.py -i input.xyz -o output.xyz -n n1 n2 n3
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
    structure = crystal()
    new_structure = crystal()
    structure.from_base_xyz(xyz)
    supercell = structure.build_supercell(args.supern)
    new_structure.get_cell_atoms(cell=supercell["cell"], atoms=supercell["atoms"])
    new_structure.to_base_xyz().to_xyz_file(args.output)
