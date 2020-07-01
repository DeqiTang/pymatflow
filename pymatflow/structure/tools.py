"""
providing tools for structure manipulation
"""

import sys
import copy
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
        normal_of_ab = np.cross(structure.cell[0], structure.cell[1])
        # get the cos of the angle between c and outer product of ab
        cosangle = np.dot(np.array(structure.cell[2]), normal_of_ab) / np.linalg.norm(np.array(structure.cell[2])) / np.linalg.norm(normal_of_ab)
        proj_c_on_ab_normal = np.linalg.norm(np.array(structure.cell[2])) * np.abs(cosangle)
        z_all = []
        for atom in structure.atoms:
            z_all.append(atom.z)
        factor = (max(z_all) - min(z_all) + thickness)  / proj_c_on_ab_normal * 1.0
        structure.cell[2] = list(np.array(structure.cell[2]) * factor)
    elif plane == 2:
        normal_of_ac = np.cross(structure.cell[0], structure.cell[2])
        # get the cos of the angle between b and outer product of ac
        cosangle = np.dot(np.array(structure.cell[1]), normal_of_ac) / np.linalg.norm(np.array(structure.cell[1])) / np.linalg.norm(normal_of_ac)
        proj_b_on_ac_normal = np.linalg.norm(np.array(structure.cell[1])) * np.abs(cosangle)
        y_all = []
        for atom in structure.atoms:
            y_all.append(atom.y)
        factor = (max(y_all) - min(y_all) + thickness)  / proj_b_on_ac_normal * 1.0
        structure.cell[1] = list(np.array(structure.cell[1]) * factor)        
    elif plane == 3:
        normal_of_bc = np.cross(structure.cell[1], structure.cell[2])
        # get the cos of the angle between a and outer product of bc
        cosangle = np.dot(np.array(structure.cell[0]), normal_of_bc) / np.linalg.norm(np.array(structure.cell[0])) / np.linalg.norm(normal_of_bc)
        proj_a_on_bc_normal = np.linalg.norm(np.array(structure.cell[0])) * np.abs(cosangle)
        x_all = []
        for atom in structure.atoms:
            x_all.append(atom.x)
        factor = (max(x_all) - min(x_all) + thickness)  / proj_a_on_bc_normal * 1.0
        structure.cell[0] = list(np.array(structure.cell[0]) * factor)
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
    
    
def enlarge_atoms(structure):
    """
    :return out:
        atoms: [
                ["C", 0.00000, 0.000000, 0.0000],
                ["O", 0.00000, 0.500000, 0.0000],
                ...
            ]
    Note: will enlarge the atoms in the unit cell along both a, b, c and -a, -b, -c direction.
        The goal is to make sure when the cell rotate in the 3D space, it will always be filled
        with atoms.
    """
    from pymatflow.base.atom import Atom
    #
    cell = copy.deepcopy(structure.cell)
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    
    n1 = np.ceil(np.max([a, b, c]) / a ) * 2 # maybe times 2 is not needed
    n2 = np.ceil(np.max([a, b, c]) / b ) * 2
    n3 = np.ceil(np.max([a, b, c]) / c ) * 2
    n = [int(n1), int(n2), int(n3)]
    print(n)
    
    atoms = copy.deepcopy(structure.atoms)
    # build supercell: replica in three vector one by one
    for i in range(3):
        natom_now = len(atoms)
        for j in range(n[i] - 1):
            for atom in atoms[:natom_now]:
                x = atom.x + float(j + 1) * structure.cell[i][0]
                y = atom.y + float(j + 1) * structure.cell[i][1]
                z = atom.z + float(j + 1) * structure.cell[i][2]
                atoms.append(Atom(atom.name, x, y, z))
        # replicate in the negative direction of structure.cell[i]
        for atom in atoms[:natom_now*n[i]]:
            x = atom.x - float(n[i]) * structure.cell[i][0]
            y = atom.y - float(n[i]) * structure.cell[i][1]
            z = atom.z - float(n[i]) * structure.cell[i][2]
            atoms.append(Atom(atom.name, x, y, z))
    return [[atom.name, atom.x, atom.y, atom.z] for atom in atoms]

def enlarge_atoms_new_cell(structure, new_cell):
    """
    :return out:
        atoms: [
                ["C", 0.00000, 0.000000, 0.0000],
                ["O", 0.00000, 0.500000, 0.0000],
                ...
            ]
    Note: will enlarge the atoms in the unit cell along both a, b, c and -a, -b, -c direction of the new_cell !!!
        but the cell is not redefined, the returned atoms is not used to form crystal, but to be 
        tailored by redefine_lattice function to get atoms for the redfined lattice.
        The goal is to make sure when the new cell rotate in the 3D space, it will always be filled
        with atoms.    
    """
    from pymatflow.base.atom import Atom
    #
    cell = copy.deepcopy(structure.cell)
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    
    new_a = np.linalg.norm(new_cell[0])
    new_b = np.linalg.norm(new_cell[1])
    new_c = np.linalg.norm(new_cell[2])
    
    n1 = np.ceil(np.max([new_a, new_b, new_c]) / a ) * 2 # maybe times 2 is not needed
    n2 = np.ceil(np.max([new_a, new_b, new_c]) / b ) * 2
    n3 = np.ceil(np.max([new_a, new_b, new_c]) / c ) * 2
    n = [int(n1), int(n2), int(n3)]
    
    atoms = copy.deepcopy(structure.atoms)
    # build supercell: replica in three vector one by one
    for i in range(3):
        natom_now = len(atoms)
        for j in range(n[i] - 1):
            for atom in atoms[:natom_now]:
                x = atom.x + float(j + 1) * structure.cell[i][0]
                y = atom.y + float(j + 1) * structure.cell[i][1]
                z = atom.z + float(j + 1) * structure.cell[i][2]
                atoms.append(Atom(atom.name, x, y, z))
        # replicate in the negative direction of structure.cell[i]
        for atom in atoms[:natom_now*n[i]]:
            x = atom.x - float(n[i]) * structure.cell[i][0]
            y = atom.y - float(n[i]) * structure.cell[i][1]
            z = atom.z - float(n[i]) * structure.cell[i][2]
            atoms.append(Atom(atom.name, x, y, z))
    return [[atom.name, atom.x, atom.y, atom.z] for atom in atoms]
    

def redefine_lattice(structure, a, b, c):
    """
    :param a, b, c: new lattice vectors in terms of old.
        new_a = a[0] * old_a + a[1] * old_b + a[2] * old_c
        like a=[1, 0, 0], b=[0, 1, 0], c=[0, 0, 1] actually defines the
        same lattice as old.
    :return an object of crystal()
    Method:
        first make a large enough supercell, which guarantee that all the atoms in the new lattice are inside
        the supercell.
        then redfine the cell, and calc the fractional coord of all atoms with regarding the new cell
        finally remove those atoms who's fractional coord is not within range [0, 1), and we can convert fractional
        coords to cartesian.
    Note:
        relationship among convertion of coords. the most important point is that all coords actually have one common
        reference system, namely the General XYZ coordinate system. all the cell are defined with XYZ as reference,
        and the convmat build from the cell(with XYZ as reference) can be applied only to atoms also with XYZ as ference,
        finally we convert frac to cartesian using convmat also defined using cell with XYZ as reference, so we get 
        the cartesian with general XYZ as reference. In the last all the coord of atoms and cell have the general 
        XYZ  system as reference. So it works!
    """
    from pymatflow.structure.crystal import crystal
    from pymatflow.base.atom import Atom
    old_cell = copy.deepcopy(structure.cell)
    new_cell = copy.deepcopy(structure.cell)
    new_cell[0] = list(a[0] * np.array(old_cell[0]) + a[1] * np.array(old_cell[1]) + a[2] * np.array(old_cell[2]))
    new_cell[1] = list(b[0] * np.array(old_cell[0]) + b[1] * np.array(old_cell[1]) + b[2] * np.array(old_cell[2]))
    new_cell[2] = list(c[0] * np.array(old_cell[0]) + c[1] * np.array(old_cell[1]) + c[2] * np.array(old_cell[2]))
    
    
    # enlarge the system
    atoms_container = crystal()
    atoms_container.get_atoms(enlarge_atoms_new_cell(structure=structure, new_cell=new_cell))
    
    # now calc the fractional coordinates of all atoms in atoms_container with new_cell as reference
    atoms_frac = []
    latcell_new = np.array(new_cell)
    convmat_new = np.linalg.inv(latcell_new.T)
    for i in range(len(atoms_container.atoms)):
        atom = []
        atom.append(atoms_container.atoms[i].name)
        atom = atom + list(convmat_new.dot(np.array([atoms_container.atoms[i].x, atoms_container.atoms[i].y, atoms_container.atoms[i].z])))
        atoms_frac.append(atom)
    
    atoms_frac_within_new_cell = []
    for atom in atoms_frac:
        if 0 <= atom[1] < 1 and 0 <= atom[2] < 1 and 0 <= atom[3] < 1:
            atoms_frac_within_new_cell.append(atom)
            
    # now convert coord of atom in atoms_frac_within_new_cell to cartesian
    out = crystal()
    out.atoms = []
    latcell_new = np.array(new_cell)
    convmat_frac_to_cartesian = latcell_new.T
    for atom in atoms_frac_within_new_cell:
        cartesian = list(convmat_frac_to_cartesian.dot(np.array([atom[1], atom[2], atom[3]])))
        out.atoms.append(Atom(name=atom[0], x=cartesian[0], y=cartesian[1], z=cartesian[2]))
    #
    
    out.cell = new_cell
    
    return out
    
    
def cleave_surface(structure, direction, thickness=10):
    """
    :param structure: an instance of crystal()
    :param direction: direction of the surface plane, like [0, 0, 1], the reference of it is three lattice 
        vector a, b, c
    
    :return an object of crystal()
    
    Note: make use of redefine_lattice() to cleave surface.
        we try to find new_a and new_b which form the surface plane (the direction is normal to the plane)
    """
    from pymatflow.structure.crystal import crystal
    #from pymatflow.base.atom import Atom
    
    old_cell = copy.deepcopy(structure.cell)
    direction_xyz_as_ref = list(direction[0] * np.array(old_cell[0]) + direction[1] * np.array(old_cell[1]) + direction[2] * np.array(old_cell[2]))
    a_from_old = None #[1, 0, 0]
    b_from_old = None #[0, 1, 0]
    c_from_old = None #[0, 0, 1]
    
    iter_ijk = []
    iter_ijk.append(0)
    for i in range(1, 16):
        iter_ijk.append(i)
        iter_ijk.append(-i)
    for i in iter_ijk:
        if a_from_old != None and b_from_old != None:
            break
        for j in iter_ijk:
            if a_from_old != None and b_from_old != None:
                break        
            for k in iter_ijk:
                if a_from_old != None and b_from_old != None:
                    break
                if i == j == k == 0:
                    continue
                new_vec = list(i * np.array(old_cell[0]) + j * np.array(old_cell[1]) + k * np.array(old_cell[2]))
                if np.dot(np.array(new_vec), np.array(direction_xyz_as_ref)) == 0:
                    if a_from_old == None:
                        a_from_old = [i, j, k]
                        new_a = list(a_from_old[0] * np.array(old_cell[0]) + a_from_old[1] * np.array(old_cell[1]) + a_from_old[2] * np.array(old_cell[2]))
                        continue
                    elif b_from_old == None:
                        cosangle = np.dot(np.array(new_vec), np.array(new_a)) / np.linalg.norm(np.array(new_vec)) / np.linalg.norm(np.array(new_a))
                        sinangle = np.sin(np.arccos(cosangle))
                        if np.abs(sinangle - 0) < 1.0e-3:
                            continue
                        elif 1.0e-5 < cosangle < 1:
                            # ignore angle between [0, 90) by default
                            # like we  ignore angle = 60 degree for a and b here, we seek for 120 instead
                            continue
                        else:
                            b_from_old = [i, j, k]
                            continue
                
    #
    new_c_ijk_sin = []
    for i in iter_ijk:
        for j in iter_ijk:
            for k in iter_ijk:
                new_c = list(i * np.array(old_cell[0]) + j * np.array(old_cell[1]) + k * np.array(old_cell[2]))
                cosangle = np.dot(np.array(new_c), np.array(direction_xyz_as_ref)) / np.linalg.norm(np.array(new_c)) / np.linalg.norm(np.array(direction_xyz_as_ref))
                sinangle = np.sin(np.arccos(cosangle))
                tmp = [i, j, k, sinangle]
                new_c_ijk_sin.append(tmp)
    sin = []
    for item in new_c_ijk_sin:
        # sometime item[3] maybe nan, we have to remove it 
        if not np.isnan(item[3]):
            sin.append(item[3])
    min_sin = min(sin)
    for item in new_c_ijk_sin:
        if item[3] == min_sin:
            c_from_old = item[:3]
            break
    print("a", a_from_old)
    print("b", b_from_old)
    print("c", c_from_old)
    # redefine lattice
    out = redefine_lattice(structure=structure, a=a_from_old, b=b_from_old, c=c_from_old)
    
    vacuum_layer(structure=out, plane=1, thickness=thickness)
    
    return out
    
def merge_layers(structure1, structure2, use_cell=None, distance=3.4, thickness=10):
    """
    :param structure1: an instance of crystal()
    :param structure2: an instance of crystal()    
    :param use_cell: use cell parameter of structure 1 or 2 or None(average) to set the new a b cell parameter
        attention: c vector is not handled this way
        
    :param distance: the distance between layers
    :param thickness: the vaccum layer thickness of the combined system
    
    :return an object of crystal()
    Note:
        only merge layers with ab plane as the surface plane
    """
    from pymatflow.structure.crystal import crystal
    from pymatflow.base.atom import Atom                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    
    old_cell_1 = copy.deepcopy(structure1.cell)
    old_cell_2 = copy.deepcopy(structure2.cell)
    
    # first transfer to fractional coordinate
    structure1.natom = len(structure1.atoms)
    frac_1 = structure1.get_fractional()
    structure2.natom = len(structure2.atoms)
    frac_2 = structure2.get_fractional()

    average_cell = []
    for i in range(3):
        average_cell.append(list((np.array(old_cell_1[i]) + np.array(old_cell_2[i])) / 2))
    
    out = crystal()
    
    if use_cell == 1:
        latcell_frac_to_cart_1 = old_cell_1
        latcell_frac_to_cart_2 = old_cell_1[0:2]
        latcell_frac_to_cart_2.append(old_cell_2[2])
    elif use_cell == 2:
        latcell_frac_to_cart_2 = old_cell_2
        latcell_frac_to_cart_1 = old_cell_2[0:2]
        latcell_frac_to_cart_1.append(old_cell_1[2])
    else:
        #average_ab = []
        #for i in range(2):
        #    vec = []
        #    for j in range(3):
        #        vec.append((old_cell_1[i][j] + old_cell_2[i][j]) / 2 )
        #    average_ab.append(vec)
        
        latcell_frac_to_cart_1 = average_cell
        latcell_frac_to_cart_1[2] = old_cell_1[2]
        latcell_frac_to_cart_2= average_cell
        latcell_frac_to_cart_2[2] = old_cell_2[2]
        
        

    # convert frac to cartesian again
    convmat_1 = np.array(latcell_frac_to_cart_1).T
    convmat_2 = np.array(latcell_frac_to_cart_2).T
    
    cart_1 = []
    for atom in frac_1:
        cartesian = list(convmat_1.dot(np.array([atom[1], atom[2], atom[3]])))
        cart_1.append([atom[0], cartesian[0], cartesian[1], cartesian[2]])
    
    cart_2 = []
    for atom in frac_2:
        cartesian = list(convmat_2.dot(np.array([atom[1], atom[2], atom[3]])))
        cart_2.append([atom[0], cartesian[0], cartesian[1], cartesian[2]])
        
    # make distance gap between cart_1 and cart_2 is the value of distance
    z_1 = []
    for atom in cart_1:
        z_1.append(atom[3])
    z_2 = []
    for atom in cart_2:
        z_2.append(atom[3])
    max_z_1 = max(z_1)
    min_z_2 = min(z_2)
    
    for i in range(len(cart_2)):
        cart_2[i][3] += distance - (min_z_2 - max_z_1)
    
    cart_all = cart_1 + cart_2
    out.atoms = []
    for atom in cart_all:
        out.atoms.append(Atom(name=atom[0], x=atom[1], y=atom[2], z=atom[3]))
    
    if use_cell == 1:
        out.cell = old_cell_1
        factor = (np.linalg.norm(np.array(old_cell_1[2])) + np.linalg.norm(np.array(old_cell_2[2]))) / np.linalg.norm(np.array(old_cell_1[2]))
        out.cell[2] = list(np.array(old_cell_1[2]) * factor)
    elif use_cell == 2:
        out.cell = old_cell_2
        factor = (np.linalg.norm(np.array(old_cell_1[2])) + np.linalg.norm(np.array(old_cell_2[2]))) / np.linalg.norm(np.array(old_cell_2[2]))
        out.cell[2] = list(np.array(old_cell_2[2]) * factor)        
    else:        
        out.cell = average_cell
        factor = (np.linalg.norm(np.array(old_cell_1[2])) + np.linalg.norm(np.array(old_cell_2[2]))) / np.linalg.norm(np.array(out.cell[2]))
        out.cell[2] = list(np.array(out.cell[2]) * factor)

    vacuum_layer(structure=out, plane=1, thickness=thickness)
    
    return out