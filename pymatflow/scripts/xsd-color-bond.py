#!/usr/bin/env python

import os
import sys
import copy
import argparse

import numpy as np
import matplotlib.pyplot as plt
from xml.etree.ElementTree import parse

from pymatflow.cmd.structflow import read_structure

"""
Color Map:
https://matplotlib.org/3.1.0/gallery/color/custom_cmap.html
"""



#from pymatflow.cmd.structflow import read_structure

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str, required=True,
        help="input structure file")

    parser.add_argument("-o", "--output", type=str, default="./contour",
        help="prefix of the output image file name")
      
    parser.add_argument("--levels", type=int, default=3, 
        help="number of levels in color map or alpha channel, etc.")
        
    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()
    # read xsd file
    xsd = parse(args.input)
    structure = read_structure(args.input)
    
    # ID of Atom3D in xsd file start from 4
    imap = xsd.getroot().find("AtomisticTreeRoot").find("SymmetrySystem").find("MappingSet").find("MappingFamily").find("IdentityMapping")
    atoms = imap.findall("Atom3d")
    bonds = imap.findall("Bond")
    natom = len(atoms)
    nbond = len(bonds)
    bonds_length = []
    for bond in bonds:
        connects = []
        connects.append(int(bond.attrib["Connects"].split(",")[0]))
        connects.append(int(bond.attrib["Connects"].split(",")[1]))
        connect_0_in_atoms = False 
        connect_1_in_atoms = False
        for i in range(len(atoms)):
            if int(atoms[i].attrib["ID"]) == connects[0]:
                connect_0_in_atoms = True
                connect_0_atom_i = i
                break
        for i in range(len(atoms)):
            if int(atoms[i].attrib["ID"]) == connects[1]:
                connect_1_in_atoms = True
                connect_1_atom_i = i
                break
        if connect_0_in_atoms == True and connect_1_in_atoms == True:            
            x0 = float(atoms[connect_0_atom_i].attrib["XYZ"].split(",")[0])
            y0 = float(atoms[connect_0_atom_i].attrib["XYZ"].split(",")[1])
            z0 = float(atoms[connect_0_atom_i].attrib["XYZ"].split(",")[2])

            x1 = float(atoms[connect_1_atom_i].attrib["XYZ"].split(",")[0])
            y1 = float(atoms[connect_1_atom_i].attrib["XYZ"].split(",")[1])
            z1 = float(atoms[connect_1_atom_i].attrib["XYZ"].split(",")[2])
           
            # transform from fractional coord to cartesian
            latcell = np.array(structure.cell)
            convmat = latcell.T
            x0, y0, z0 = list(convmat.dot(np.array([x0, y0, z0])))        
            x1, y1, z1 = list(convmat.dot(np.array([x1, y1, z1])))
            length = np.sqrt((x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2)
            bonds_length.append(length)
        else:
            bonds_length.append(None)
    #   
    bonds_length_norm = []
    tmp = copy.deepcopy(bonds_length)
    while None in tmp:
        tmp.remove(None)
    #print(tmp)
    bond_max = max(tmp)
    bond_min = min(tmp)
    print(bond_max, bond_min)
    for i in range(len(bonds_length)):
        if bonds_length[i] == None:
            bonds_length_norm.append(None)
        else:
            bonds_length_norm.append((bonds_length[i] - bond_min) / (bond_max - bond_min))
            
    print(bonds_length_norm)
    # ------------------------------------------------------------------------
    # now bonds length are alreay handled
    # we should use different scheme to set the color of them
    # ------------------------------------------------------------------------
    
    #Scheme I
    # set vriridis cmap color
    from matplotlib import cm
    #cmap = cm.get_cmap("viridis", args.levels)
    #cmap = cm.get_cmap("Blues", args.levels)
    #cmap = cm.get_cmap("Greens", args.levels)
    #cmap = cm.get_cmap("Reds", args.levels)
    #cmap = cm.get_cmap("gist_heat", args.levels)
    #cmap = cm.get_cmap("hsv", args.levels)
    cmap = cm.get_cmap("rainbow", args.levels)
    #cmap = cm.get_cmap("winter", args.levels)
    for i in range(len(bonds)):
        if bonds_length_norm[i] == None:
            pass
        else:
            # there are two ways to get the color for one value, if you use n color in cmap
            # if you pass a integer (0<=value<=(n-1)) to cmap(), you will get the color of it
            # or you can pass a float within([0.0, 1.0]) to cmap(), you will get the color of it
            # as bonds_length_norm is normalized to (0.0-1.0), we could directly pass it to cmap() to get the right color
            color = cmap(bonds_length_norm[i])
            # cmap return RGB in range of (0, 1), we have to multiply it with 255 to get a value between range(0, 255)
            bonds[i].set("Color", "%f, %f, %f, %f" % (color[0]*255, color[1]*255, color[2]*255, color[3]))
    # plot the cmap
    # the plot is with no use, what we want is the corresponding color bar.
    tmp_norm = copy.deepcopy(bonds_length_norm)
    while None in tmp_norm:
        tmp_norm.remove(None)
    tmp_value = [tmp_norm[i] * (bond_max - bond_min) + bond_min for i in range(len(tmp_norm))]
    color_bar = plt.scatter(tmp_norm, tmp_norm, c=tmp_value, cmap=cmap)
    plt.colorbar(color_bar)
    plt.savefig(args.output+".colorbar.png")
    plt.close()
    
        
    """
    #Scheme II
    for i in range(len(bonds)):
        if bonds_length_norm[i] == None:
            pass
        else:
            bonds[i].set("Color", "%f, %f, %f, %f" % (35.3, 35.3, 35.3, bonds_length_norm[3]))    
    """
    
    
    
    
    # write xsd file
    xsd.write(args.output)



if __name__ == "__main__":
    main()