#!/usr/bin/env python

import os
import argparse
import numpy as np

from pymatflow.structure.crystal import crystal
"""
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--xdatcar", type=str,
            help="the XDATCAR file")

    parser.add_argument("-o", "--xyz", type=str,
            help="the output xyz file")

    args = parser.parse_args()

    print("========================================\n")
    print("         XDATCAR-to-xyz.py\n")
    print("----------------------------------------\n")

    with open(args.xdatcar, 'r') as fin:
        xdatcar = fin.readlines()
        
    natom = 0
    for i in xdatcar[6].split("\n")[0].split():
        natom += int(i)
        
    nimage = 0
    for line in xdatcar:
        if line.split()[0] == "Direct":
            nimage += 1
    
    direct_begin_lines = []
    for i in range(len(xdatcar)):
        if xdatcar[i].split()[0] == "Direct":
            direct_begin_lines.append(i)
    
        
    if nimage < 2:
        vc = False
    else:
        if direct_begin_lines[1] - direct_begin_lines[0] != natom + 1:
            vc = True
        else:
            vc = False
    
    element_n_list = []
    for i in range(len(xdatcar[5].split())):
        element = {}
        element["label"] = xdatcar[5].split()[i]
        element["n"] = int(xdatcar[6].split()[i])
        element_n_list.append(element)
    
    element_list = []
    for i in range(len(element_n_list)):
        for j in range(element_n_list[i]["n"]):
            element_list.append(element_n_list[i]["label"])

    
    images = []
    
    for i in range(nimage):
        image = crystal()
        image.cell = []
        if i == 0:
            for j in range(3):
                tmp = [float(xdatcar[2+j].split()[0]), float(xdatcar[2+j].split()[1]), float(xdatcar[2+j].split()[2])]
                image.cell.append(tmp)
        else:
            if vc == False:
                image.cell = images[-1].cell
            else:
                for j in range(3):
                    tmp = [float(xdatcar[direct_begin_lines[i]-5+j].split()[0]), float(xdatcar[direct_begin_lines[i]-5+j].split()[1]), float(xdatcar[direct_begin_lines[i]-5+j].split()[2])]
                    image.cell.append(tmp)
        atoms_list = []
        latcell = np.array(image.cell)
        convmat = latcell.T
        for j in range(natom):
            cartesian = list(convmat.dot(np.array([
                float(xdatcar[direct_begin_lines[i]+j+1].split()[0]), 
                float(xdatcar[direct_begin_lines[i]+j+1].split()[1]), 
                float(xdatcar[direct_begin_lines[i]+j+1].split()[2])
            ])))
        
            atoms_list.append([
                element_list[j],
                cartesian[0],
                cartesian[1],
                cartesian[2]
            ])
        image.get_atoms(atoms_list)
        images.append(image)
        
    # write structure file xyz
    with open(args.xyz, 'w') as fout:
        for image in images:
            fout.write("%d\n" % len(image.atoms))
            fout.write("cell: %f %f %f | %f %f %f | %f %f %f\n" % (
                image.cell[0][0], image.cell[0][1], image.cell[0][2],
                image.cell[1][0], image.cell[1][1], image.cell[1][2],
                image.cell[2][0], image.cell[2][1], image.cell[2][2],
            ))
            for atom in image.atoms:
                fout.write("%s\t%f\t%f\t%f\n" % (
                    atom.name,
                    atom.x,
                    atom.y,
                    atom.z
                ))
    # end write structure file xyz