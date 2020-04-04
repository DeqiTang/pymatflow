"""
providing tools for structure manipulation
"""

import sys
import numpy as np

def move_along(structure, atoms_to_move, direc, disp):
    """
    :param structure: an instance of pymatflow.structure.crystal.crystal()
    :params atoms_to_move: a list of atoms to move, counting starts with 0
    :param direc: three number to indicate the direction to move along
        namely the crystal orientation index
    :param disp: displacement of the atoms in unit of Angstrom
    """
    # do some checking
    if max(atoms_to_move) > (len(structure.atoms) - 1):
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

    for i in atoms_to_move:
        structure.atoms[i].x += deltax
        structure.atoms[i].y += deltay
        structure.atoms[i].z += deltaz
    # end

def remove_atoms(structure, atoms_to_remove):
    """
    :param structure: an instance of pymatflow.structure.crystal.crystal()
    :params atoms_to_remove: a list of atoms to remove, counting starts with 0
    """
    # do some checking
    if max(atoms_to_remove) > (len(structure.atoms) - 1):
        print("=============================================================================\n")
        print("                              WARNING \n")
        print("------------------------------------------------------------------------------\n")
        print("the atom you are trying to remove is beyond the number of atoms in the structure\n")
        sys.exit(1)

    for i in atoms_to_remove:
        structure.atoms[i] = None
    while None in structure.atoms:
        structure.atoms.remove(None)
    # end    

def vacuum_layer(structure, plane, thickness):
    """
    :param structure: an instance of pymatflow.structure.crystal.crystal()
    :param plane: the plane chosen to add vacuum layer, can be:
        1 -> ab plane
        2 -> ac plane
        3 -> bc plane
    :param thickness: thickness of the vacuum layer
    """
    if plane == 1:
        ab = np.cross(structure.cell[0], structure.cell[1])
        # get the cos of the angle between c and outer product of ab
        cos_angle = np.dot(structure.cell[2], ab) / (np.linalg.norm(structure.cell[2]), np.linalg.norm(ab))
        scale_c = (thickness + np.linalg.norm(structure.cell[2])) / np.linalg.norm(structure.cell[2])
        for i in range(3):
            structure.cell[2][i] *= scale_c
    elif plane == 2:
        ac = np.cross(structure.cell[0], structure.cell[2])
        # get the cos of the angle between c and outer product of ab
        cos_angle = np.dot(structure.cell[1], ac) / (np.linalg.norm(structure.cell[1]), np.linalg.norm(ac))
        scale_b = (thickness + np.linalg.norm(structure.cell[1])) / np.linalg.norm(structure.cell[1])
        for i in range(3):
            structure.cell[1][i] *= scale_b            
    elif plane == 3:
        bc = np.cross(structure.cell[1], structure.cell[2])
        # get the cos of the angle between c and outer product of ab
        cos_angle = np.dot(structure.cell[0], bc) / (np.linalg.norm(structure.cell[0]), np.linalg.norm(bc))
        scale_a = (thickness + np.linalg.norm(structure.cell[0])) / np.linalg.norm(structure.cell[0])
        for i in range(3):
            structure.cell[0][i] *= scale_a
    else:
        pass
    # end