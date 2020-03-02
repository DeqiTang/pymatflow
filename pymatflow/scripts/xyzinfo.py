#!/usr/bin/env python

import numpy as np
from pymatflow.base.xyz import base_xyz
import seekpath
import spglib
import argparse

"""
Reference:
    https://atztogo.github.io/spglib/python-spglib.html
    https://seekpath.readthedocs.io/en/latest/index.html
Warning:
    the result is not guaranteed to be correct
"""

def get_spacegroup_and_kpath(xyz):
    """
    xyz is an instance of base_xyz
    """
    lattice = xyz.cell  # [xyz.cell[0:3], xyz.cell[3:6], xyz.cell[6:9]]
    positions = []
    numbers = []
    #a = np.sqrt(xyz.cell[0]**2 + xyz.cell[1]**2 + xyz.cell[2]**2)
    #b = np.sqrt(xyz.cell[3]**2 + xyz.cell[4]**2 + xyz.cell[5]**2)
    #c = np.sqrt(xyz.cell[6]**2 + xyz.cell[7]**2 + xyz.cell[8]**2)
    a = np.sqrt(xyz.cell[0][0]**2 + xyz.cell[0][1]**2 + xyz.cell[0][2]**2)
    b = np.sqrt(xyz.cell[1][0]**2 + xyz.cell[1][1]**2 + xyz.cell[1][2]**2)
    c = np.sqrt(xyz.cell[2][0]**2 + xyz.cell[2][1]**2 + xyz.cell[2][2]**2)
    for atom in xyz.atoms:
        positions.append([atom.x / a, atom.y / b, atom.z / c]) # must be scaled cartesian
        numbers.append(xyz.specie_labels[atom.name])
    cell = (lattice, positions, numbers)
    return [spglib.get_spacegroup(cell, symprec=1.0e-5), seekpath.get_path(cell)]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str,
            help="the xyz structure file")
    args = parser.parse_args()


    xyz = base_xyz()
    xyz.get_xyz(args.file)

    print("===========================================\n")
    print("calculated using spglib\n")
    print("===========================================\n")
    print("spacegroup is : %s\n" % get_spacegroup_and_kpath(xyz)[0])
    print("Warning:\n")
    print("the result is not guaranteed to be correct\n")
    print("-------------------------------------------\n")
    print("suggested k path calculated using seekpath:\n")
    print(get_spacegroup_and_kpath(xyz)[1])
