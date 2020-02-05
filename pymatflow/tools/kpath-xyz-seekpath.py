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
"""

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--xyz", type=str, default=None,
            help="the input xyz structure file")
    
    parser.add_argument("-o", "--output", type=str, default="kpath-from-seekpath.txt",
            help="the output kpoitns file")

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
    
    specialk = [] 
    # [[kx, ky, kz, label, end_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', None], ...]
    # if end_indicator in a qpoint is None, then it will connect to the following point,
    # if end_indicator in a qpoint is '|', then it will not connect to the following point,
    specialk.append([
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][0]][0]),
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][0]][1]),
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][0]][2]),
        kpoints_seekpath["path"][0][0],
        None,
        ])
    specialk.append([
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][1]][0]),
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][1]][1]),
        float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][1]][2]),
        kpoints_seekpath["path"][0][1],
        None, # first set to None here, and it will be reset by the following path by judge whether it is connected 
        ])
    for i in range(1, len(kpoints_seekpath["path"])):
        if kpoints_seekpath["path"][i][0] == kpoints_seekpath["path"][i-1][1]:
            specialk[-1][4] = None
            specialk.append([
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][0]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][1]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][2]),
                kpoints_seekpath["path"][i][1],
                None,
                ])
        else:
            specialk[-1][4] = "|"
            specialk.append([
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][0]][0]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][0]][1]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][0]][2]),
                kpoints_seekpath["path"][i][0],
                None,
                ])
            specialk.append([
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][0]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][1]),
                float(kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]][2]),
                kpoints_seekpath["path"][i][1],
                None,
                ])
    #
    with open(args.output, 'w') as fout:
        fout.write("%d\n" % len(specialk))
        for kpoint in specialk:
            if kpoint[4] == None:
                fout.write("%f %f %f #%s\n" % (kpoint[0], kpoint[1], kpoint[2], kpoint[3]))
            elif kpoint[4] == "|":
                fout.write("%f %f %f #%s |\n" % (kpoint[0], kpoint[1], kpoint[2], kpoint[3]))
    #
    

    print("===========================================\n")
    print("calculated using seekpath\n")
    print("===========================================\n")
    print("Warning:\n")
    print("the result is not guaranteed to be correct\n")
    print("-------------------------------------------\n")
    print("suggested k path calculated using seekpath:\n")
    print(kpoints_seekpath)
