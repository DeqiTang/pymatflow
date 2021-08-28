import copy
import numpy as np
import seekpath
import spglib

from pymatflow.base.element import element
"""
Reference:
    https://atztogo.github.io/spglib/python-spglib.html
    https://seekpath.readthedocs.io/en/latest/index.html
Warning:
    the result is not guaranteed to be correct

    seekpath will automatically find the primitive cell of the structure
    input by you. and the k points it gives is corresponding with that
    primitive cell, if you stick to the use the original structure, you
    must know how to modify the k points to be applicable to your original
    structure.

    seekpath generated high symmetry are in crystal coordinates not cartesian.
    so we should set K_POINTS in qe to {crystal_b}
    p.s. like there are cartesian and crystal(fractional) coordinates in real
    space for structure coordinates, there are also cartesian and crystal coordinates
    for reciprocal space for kpoint. in qe, tpiba and tpiba_b are cartesian coordinates
    in unit of 2pi/a, the latter with suffix '_b' means 'for band structure'.
    and crystal and crysta_b are in crystal coordinated. seek-path generated k points
    are in reciprocal crystal coordinate, so we should use crystal_b(qe) or B_VECTOR(CP2K)
    for band structure calculation using seekpath generated high symmetry kpoing.
"""

def seekpath_output_kpath(structure, output, with_time_reversal=True, symprec=1.0e-9, join=15):
    """
    structure is an instance of BaseXyz or Crystal
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

    kpoints_seekpath = seekpath.get_path(
        structure=cell,
        with_time_reversal=with_time_reversal,
        symprec = symprec
        )

    kpath = []
    # [[kx, ky, kz, label, connect_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', 15], ...]
    # if connect_indicator in a kpoint is an integer, then it will connect to the following point
    # through the number of kpoints defined by connect_indicator.
    # if connect_indicator in a kpoint is '|', then it will not connect to the following point,
    kpath.append([
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][0]][0]),
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][0]][1]),
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][0]][2]),
        kpoints_seekpath["path"][0][0],
        join,
        ])
    kpath.append([
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][1]][0]),
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][1]][1]),
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][1]][2]),
        kpoints_seekpath["path"][0][1],
        join, # first set to oin here, and it will be reset by the following path by judge whether it is connected
        ])
    for i in range(1, len(kpoints_seekpath["path"])):
        if kpoints_seekpath["path"][i][0] == kpoints_seekpath["path"][i-1][1]:
            kpath[-1][4] = join
            kpath.append([
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][0]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][1]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][2]),
                kpoints_seekpath["path"][i][1],
                join,
                ])
        else:
            kpath[-1][4] = "|"
            kpath.append([
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][0]][0]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][0]][1]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][0]][2]),
                kpoints_seekpath["path"][i][0],
                join,
                ])
            kpath.append([
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][0]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][1]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][2]),
                kpoints_seekpath["path"][i][1],
                join,
                ])
    #
    with open(output, 'w') as fout:
        fout.write("%d\n" % len(kpath))
        for kpoint in kpath:
            fout.write("%f %f %f #%s %s\n" % (kpoint[0], kpoint[1], kpoint[2], kpoint[3], str(kpoint[4])))
    #
    print("===========================================\n")
    print("calculated using seekpath\n")
    print("===========================================\n")
    print("Warning:\n")
    print("the result is not guaranteed to be correct\n")
    print("-------------------------------------------\n")
    print("suggested k path calculated using seekpath:\n")
    print(kpoints_seekpath)


def get_spacegroup_and_kpath(structure, with_time_reversal=True, symprec=1.0e-9):
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

    return [
        spglib.get_spacegroup(cell, symprec=symprec),
        seekpath.get_path(structure=cell, with_time_reversal=with_time_reversal, symprec=symprec
    )]


def seekpath_std_structure(structure, with_time_reversal=True, symprec=1.0e-9):
    """
    :param: structure is an instance of crystal()
    
    Note: 
        standardize the structure using seekpath
    """
    from pymatflow.base.atom import Atom
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
    
    result = seekpath.get_path(structure=cell, with_time_reversal=with_time_reversal, symprec=symprec)

    out = copy.deepcopy(structure)
    out.cell = result["conv_lattice"]

    std_fractional = []
    for i in range(len(structure.atoms)):
        std_fractional.append(result["conv_positions"][i])
    std_cartesian = []
    
    out.atoms = []
    latcell = np.array(out.cell)
    convmat_frac_to_cartesian = latcell.T
    for i, frac in enumerate(std_fractional):
        cartesian = list(convmat_frac_to_cartesian.dot(np.array([frac[0], frac[1], frac[2]])))
        out.atoms.append(Atom(name=structure.atoms[i].name, x=cartesian[0], y=cartesian[1], z=cartesian[2]))
    #
    print(result)
    #
    return out