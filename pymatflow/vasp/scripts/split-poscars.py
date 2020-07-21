#!/usr/bin/env python

import os
import argparse
import numpy as np

from pymatflow.structure.crystal import crystal
"""
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--poscars", type=str, required=True,
            help="the file containing multi POSCAR structure")

    parser.add_argument("-d", "--directory", type=str,
            default="poscar_structures_split",
            help="directory to put the splitted structures.")

    args = parser.parse_args()

    print("========================================\n")
    print("         split-poscars.py\n")
    print("----------------------------------------\n")

    with open(args.poscars, 'r') as fin:
        poscars = fin.readlines()
        
    natom = 0
    for i in poscars[6].split("\n")[0].split():
        natom += int(i)
        
    nimage = 0
    for line in poscars:
        if line.split()[0] == "Direct":
            nimage += 1
    
    direct_begin_lines = []
    for i in range(len(poscars)):
        if poscars[i].split()[0] == "Direct":
            direct_begin_lines.append(i)
    
        
    if nimage < 2:
        vc = False
    else:
        if direct_begin_lines[1] - direct_begin_lines[0] != natom + 1:
            vc = True
        else:
            vc = False
    
    element_n_list = []
    for i in range(len(poscars[5].split())):
        element = {}
        element["label"] = poscars[5].split()[i]
        element["n"] = int(poscars[6].split()[i])
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
                tmp = [float(poscars[2+j].split()[0]), float(poscars[2+j].split()[1]), float(poscars[2+j].split()[2])]
                image.cell.append(tmp)
        else:
            if vc == False:
                image.cell = images[-1].cell
            else:
                for j in range(3):
                    tmp = [float(poscars[direct_begin_lines[i]-5+j].split()[0]), float(poscars[direct_begin_lines[i]-5+j].split()[1]), float(poscars[direct_begin_lines[i]-5+j].split()[2])]
                    image.cell.append(tmp)
        atoms_list = []
        latcell = np.array(image.cell)
        convmat = latcell.T
        for j in range(natom):
            cartesian = list(convmat.dot(np.array([
                float(poscars[direct_begin_lines[i]+j+1].split()[0]), 
                float(poscars[direct_begin_lines[i]+j+1].split()[1]), 
                float(poscars[direct_begin_lines[i]+j+1].split()[2])
            ])))
        
            atoms_list.append([
                element_list[j],
                cartesian[0],
                cartesian[1],
                cartesian[2]
            ])
        image.get_atoms(atoms_list)
        images.append(image)
        
    # write structure files
    os.system("mkdir -p %s" % args.directory)
    from pymatflow.cmd.structflow import write_structure
    for i in range(len(images)):
        write_structure(structure=images[i], filepath=os.path.join(args.directory, "%d.cif" % (i+1)))
        
    # end write structure files