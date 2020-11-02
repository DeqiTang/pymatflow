#!/usr/bin/env python

import os
import sys
import copy
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

from pymatflow.cmd.structflow import read_structure
from pymatflow.cmd.structflow import write_structure

from pymatflow.vasp.post.pdos import post_pdos

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str, required=True,
        help="input vasp *CHG* file, -i *CHG* ")

    parser.add_argument("--output-structure", type=str, default="chg",
        help="output stucture contained in *CHG*")

    parser.add_argument("-o", "--output", type=str, default="chg",
        help="prefix of the output image file name")
    
    parser.add_argument("--levels", type=int, default=10,
        help="levels of the color map or color bar")

    parser.add_argument("-z", "--z", type=float, default=1,
        help="a value between 0 and 1, indicat height of in z direction to print the plot")
        
    parser.add_argument("--cmap", type=str, default="gray",
        choices=["gray", "hot", "afmhot", "Spectral", "plasma", "magma", "hsv", "rainbow", "brg"])
        
    parser.add_argument("--abscissa", type=str, nargs="+", default=["a", "b", "c"], 
        choices=["a", "b", "c"], 
        help="choose the direction to do the dimension reduction")
        
    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()
    
    chg_filepath = args.input
    
    with open(chg_filepath, "r") as fin:
        chg = fin.readlines()
    
    for j in range(len(chg)):
        if len(chg[j].split()) == 0:
            first_blank_line = j
            break
    
    first_augmentation_line = None
    for j in range(len(chg)):
        if "augmentation" in chg[j]:
            first_augmentation_line = j
            break
                
    os.system("mkdir -p /tmp/pymatflow/")
    with open("/tmp/pymatflow/POSCAR", "w") as fout:
        for j in range(first_blank_line):
            fout.write(chg[j])
            
    structure = read_structure("/tmp/pymatflow/POSCAR")
    write_structure(structure=structure, filepath=args.output_structure+".cif")

    a = np.linalg.norm(structure.cell[0])
    b = np.linalg.norm(structure.cell[1])
    c = np.linalg.norm(structure.cell[2])
    
    # assume three *CHG* have the same ngxf and ngyf ngzf
    ngxf = int(chg[first_blank_line+1].split()[0])
    ngyf = int(chg[first_blank_line+1].split()[1])
    ngzf = int(chg[first_blank_line+1].split()[2])    
    

    if first_augmentation_line == None:
        tmp_str = "".join(chg[first_blank_line+2:])
        data = np.fromstring(tmp_str, sep="\n").reshape(ngzf, ngyf, ngxf)
    #data = data.reshape(ngzf, ngyf, ngxf)

    
    # data dimension reduction
    # the unit of value is actually not physical now!
    cell_volume = np.dot(np.cross(np.array(structure.cell[0]), np.array(structure.cell[1])), np.array(structure.cell[2]))
    cell_volume_per_unit = cell_volume / (ngzf * ngyf * ngxf)
    
    # value in Vasp *CHG* are \rho(r)_of_electrons * Volume_of_cell, so we should divide it by cell_volume here and time it with cell_volume_per_unit
    # to get the number of electrons per divided unit
    total_electrons = np.sum(data) / cell_volume * cell_volume_per_unit
    #
    
    print("======================================================\n")
    print("           Information collected\n")
    print("------------------------------------------------------\n")
    print("cell volume: %f (A^3)\n" % cell_volume)
    print("total electrons: %f\n" % total_electrons)
    
    
    
    # unit of data_red_? is e/Anstrom, namely number of electrons per Angstrom
    data_red_a = []
    data_red_b = []
    data_red_c = []
    if "c" in args.abscissa:
        factor = cell_volume_per_unit / cell_volume
        len_ci = c / ngzf
        for ci in range(data.shape[0]):
            tmp = 0
            for bi in range(data.shape[1]):
                tmp += np.sum(data[ci, bi, :])
            nelect_ci = tmp * factor
            rho_line = nelect_ci / len_ci
            data_red_c.append(rho_line)
    if "b" in args.abscissa:
        factor = cell_volume_per_unit / cell_volume
        len_bi = b / ngyf    
        for bi in range(data.shape[1]):
            tmp = 0
            for ai in range(data.shape[2]):
                tmp += np.sum(data[:, bi, ai])
            nelect_bi = tmp * factor
            rho_line = nelect_bi / len_bi                
            data_red_b.append(rho_line)
    if "a" in args.abscissa:
        factor = cell_volume_per_unit / cell_volume
        len_ai = a / ngxf        
        for ai in range(data.shape[2]):
            tmp = 0
            for ci in range(data.shape[0]):
                tmp += np.sum(data[ci, :, ai])
            nelect_ai = tmp * factor
            rho_line = nelect_ai / len_ai                   
            data_red_a.append(rho_line)    

    # output the data and make the plot
    if "c" in args.abscissa:
        with open(args.output+".1d.c.data", 'w') as fout:
            fout.write("#c(angstrom) rho(e) (number of electron per Angstrom)\n")
            c_coord = np.linspace(0, c, len(data_red_c))
            for i in range(len(data_red_c)):
                fout.write("%f %f\n" % (c_coord[i], data_red_c[i]))
        plt.plot(np.linspace(0, c, len(data_red_c)), data_red_c)                
        plt.ylabel(r"$\rho (e/\AA)$")
        plt.tight_layout()
        plt.savefig(args.output+".1d.c.png")
        plt.close()                
    if "b" in args.abscissa:
        with open(args.output+".1d.b.data", 'w') as fout:
            fout.write("#b(angstrom) rho(e) (number of electron per Angstrom)\n")
            b_coord = np.linspace(0, b, len(data_red_b))
            for i in range(len(data_red_b)):
                fout.write("%f %f\n" % (b_coord[i], data_red_b[i]))        
        plt.plot(np.linspace(0, b, len(data_red_b)), data_red_b)    
        plt.ylabel(r"$\rho (e/\AA)$")              
        plt.tight_layout()
        plt.savefig(args.output+".1d.b.png")
        plt.close()                
    if "a" in args.abscissa:
        with open(args.output+".1d.a.data", 'w') as fout:
            fout.write("#a(angstrom) rho(e) (number of electron per Angstrom)\n")
            a_coord = np.linspace(0, a, len(data_red_a))
            for i in range(len(data_red_a)):
                fout.write("%f %f\n" % (a_coord[i], data_red_a[i]))
        plt.plot(np.linspace(0, a, len(data_red_a)), data_red_a)                
        plt.ylabel(r"$\rho (e/\AA)$")            
        plt.tight_layout()
        plt.savefig(args.output+".1d.a.png")
        plt.close()
    

if __name__ == "__main__":
    main()