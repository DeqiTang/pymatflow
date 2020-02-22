#!/usr/bin/env python

import numpy as np
import seekpath
import argparse

from pymatflow.base.xyz import base_xyz

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

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--xyz", type=str, default=None,
            help="the input xyz structure file")
    
    parser.add_argument("-o", "--output", type=str, default="kpath-from-seekpath.txt",
            help="the output kpoitns file")

    parser.add_argument("--join", type=int, default=15,
            help="default number of kpoint to connect the connected high symmetry k point")

    args = parser.parse_args()
    #

    xyz = base_xyz()
    xyz.get_xyz(args.xyz)

    lattice = xyz.cell
    positions = []
    numbers = []
    a = np.sqrt(xyz.cell[0][0]**2 + xyz.cell[0][1]**2 + xyz.cell[0][2]**2)
    b = np.sqrt(xyz.cell[1][0]**2 + xyz.cell[1][1]**2 + xyz.cell[1][2]**2)
    c = np.sqrt(xyz.cell[2][0]**2 + xyz.cell[2][1]**2 + xyz.cell[2][2]**2)
    for atom in xyz.atoms:
        positions.append([atom.x / a, atom.y / b, atom.z / c]) # must be scaled cartesian
        numbers.append(xyz.specie_labels[atom.name])
    cell = (lattice, positions, numbers)

    kpoints_seekpath = seekpath.get_path(cell)
    
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
        args.join,
        ])
    kpath.append([
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][1]][0]),
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][1]][1]),
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][1]][2]),
        kpoints_seekpath["path"][0][1],
        args.join, # first set to args.join here, and it will be reset by the following path by judge whether it is connected 
        ])
    for i in range(1, len(kpoints_seekpath["path"])):
        if kpoints_seekpath["path"][i][0] == kpoints_seekpath["path"][i-1][1]:
            kpath[-1][4] = args.join 
            kpath.append([
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][0]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][1]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][2]),
                kpoints_seekpath["path"][i][1],
                args.join,
                ])
        else:
            kpath[-1][4] = "|"
            kpath.append([
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][0]][0]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][0]][1]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][0]][2]),
                kpoints_seekpath["path"][i][0],
                args.join,
                ])
            kpath.append([
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][0]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][1]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][2]),
                kpoints_seekpath["path"][i][1],
                args.join,
                ])
    #
    with open(args.output, 'w') as fout:
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
