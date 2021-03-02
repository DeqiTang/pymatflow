#!/usr/bin/env python

import numpy as np
from pymatflow.base.xyz import base_xyz
from pymatflow.cmd.structflow import read_structure
from pymatflow.base.element import element
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

def get_spacegroup_and_kpath(structure, symprec=1.0e-9):
    """
    :param: structure is an instance of crystal()
    """
    lattice = structure.cell
    positions = []
    numbers = []
    a = np.sqrt(structure.cell[0][0]**2 + structure.cell[0][1]**2 + structure.cell[0][2]**2)
    b = np.sqrt(structure.cell[1][0]**2 + structure.cell[1][1]**2 + structure.cell[1][2]**2)
    c = np.sqrt(structure.cell[2][0]**2 + structure.cell[2][1]**2 + structure.cell[2][2]**2)
        
    frac_coord = structure.get_fractional()
    for i in range(len(frac_coord)):
        positions.append([frac_coord[i][1], frac_coord[i][2], frac_coord[i][3]]) # must be fractional coordinates
        numbers.append(element[frac_coord[i][0]].number)

    cell = (lattice, positions, numbers)

    return [spglib.get_spacegroup(cell, symprec=symprec), seekpath.get_path(cell)]

    
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str,
            help="the input structure file")
            
    parser.add_argument("--symprec", type=float, default=1.0e-7,
            help="symmetry precision")
            
    args = parser.parse_args()

    a = read_structure(filepath=args.input)

    print("===========================================\n")
    print("calculated using spglib\n")
    print("===========================================\n")
    print("spacegroup is : %s\n" % get_spacegroup_and_kpath(structure=a, symprec=args.symprec)[0])
    print("Warning:\n")
    print("the result is not guaranteed to be correct\n")
    print("-------------------------------------------\n")
    print("suggested k path calculated using seekpath:\n")
    print(get_spacegroup_and_kpath(structure=a, symprec=args.symprec)[1])
