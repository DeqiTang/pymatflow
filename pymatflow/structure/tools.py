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

    
def inverse_geo_center(structure):
    """
    calc the geometric center of the system and make an inversion against that center
    :param structure: an instance of pymatflow.structure.crystal.crystal()
    """
    # calc the geometric center
    x = 0
    y = 0
    z = 0
    for atom in structure.atoms:
        x += atom.x
        y += atom.y
        z += atom.z
    x /= len(structure.atoms)
    y /= len(structure.atoms)
    z /= len(structure.atoms)
    # now get the symmetry image against the geometric center
    for atom in structure.atoms:
        atom.x = x * 2 - atom.x
        atom.y = y * 2 - atom.y
        atom.z = z * 2 - atom.z
    # end

def inverse_point(structure, point):
    """
    calc the geometric center of the system and make an inversion against that center
    :param structure: an instance of pymatflow.structure.crystal.crystal()
    :param point: the inverse center point, like [0.0, 0.0, 0.0]
    """
    # now get the symmetry image against the inverse center
    for atom in structure.atoms:
        atom.x = point[0] * 2 - atom.x
        atom.y = point[1] * 2 - atom.y
        atom.z = point[2] * 2 - atom.z
    # end

def inverse_cell_center(structure):
    """
    make an inversion against the cell center
    :param structure: an instance of pymatflow.structure.crystal.crystal()
    """
    # first transfer to fractional coordinate and inverse against [0.5, 0.5, 0.5]
    structure.natom = len(structure.atoms)
    frac = structure.get_fractional()
    for atom in frac:
        atom[1] = 0.5 * 2 - atom[1]
        atom[2] = 0.5 * 2 - atom[2]
        atom[3] = 0.5 * 2 - atom[3]
    # convert frac to cartesian again
    latcell = np.array(structure.cell)
    convmat = latcell.T
    from pymatflow.base.atom import Atom
    structure.atoms = []
    for atom in frac:
        cartesian = list(convmat.dot(np.array([atom[1], atom[2], atom[3]])))
        structure.atoms.append(Atom(name=atom[0], x=cartesian[0], y=cartesian[1], z=cartesian[2]))
    #

def rotate_along_axis(structure, rotate_atoms=[], axis=[]):
    """
    rotate the specified atoms along the specified axis
    :param structure: an instance of pymatflow.structure.crystal.crystal()
    """
    pass