"""
providing tools for structure manipulation
"""

import sys
import numpy as np

def move_along(structure, atom_to_move, direc, disp):
    """
    :param structure: an instance of pymatflow.structure.crystal.crystal()
    :params atom_to_move: a list of atoms to move, counting starts with 0
    :param direc: three number to indicate the direction to move along
        namely the crystal orientation index
    :param disp: displacement of the atoms in unit of Angstrom
    """
    # do some checking
    if max(atom_to_move) > (len(structure.atoms) - 1):
        print("=============================================================================\n")
        print("                              WARNING \n")
        print("------------------------------------------------------------------------------\n")
        print("the atom you are trying to move is beyond the number of atoms in the structure\n")
        sys.exit(1)

    direc_cartesian = np.array(structure.cell[0]) * direc[0] + np.array(structure.cell[1]) * direc[1] + np.array(structure.cell[2]) * direc[2]
    # normalize
    length = np.sqrt(direc_cartesian[0]**2+direc_cartesian[1]**2+direc_cartesian[2]**2)
    direc_cartesian = direc_cartesian / length

    deltax = direc_cartesian[0] * disp
    deltay = direc_cartesian[1] * disp
    deltaz = direc_cartesian[2] * disp

    for i in atom_to_move:
        structure.atoms[i].x += deltax
        structure.atoms[i].y += deltay
        structure.atoms[i].z += deltaz
    # end