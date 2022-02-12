#!/usr/bin/env python

import os
import sys
import copy
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms


from pymatflow.base.element import element as elem
from pymatflow.structure.crystal import crystal


from pymatflow.cmd.structflow import read_structure
from pymatflow.cmd.structflow import write_structure


def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str, nargs=3, required=True,
        help="input  cube file, -i TOTAL PART1 PART2")

    parser.add_argument("--output-structure", type=str, default="diff-chg",
        help="output stucture contained in cube file")

    parser.add_argument("-o", "--output", type=str, default="diff-chg",
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
    
    cube_filepath = args.input
    
    cube = []
    for i in range(3):
        with open(cube_filepath[i], "r") as fin:
            cube.append(fin.readlines())

    bohr_to_angstrom = 0.529177249
    natoms_each_image = []
    for k in range(3):
        # read structure info
        natom = abs(int(cube[k][2].split()[0])) # it might be negative, if MO infor are included in cube file
        natoms_each_image.append(natom)
        structure = crystal()
        structure.cell = []
        for i in range(3):
            tmp = []
            for j in range(3):
                tmp.append(int(cube[k][i+3].split()[0]) * float(cube[k][i+3].split()[j+1]) * bohr_to_angstrom)
            structure.cell.append(tmp)
        atoms_list = []
        for i in range(natom):
            atomic_number = int(cube[k][i+6].split()[0])
            for e in elem:
                if elem[e].number == atomic_number:
                    label = e
            atoms_list.append([
                label,
                float(cube[k][i+6].split()[2]) * bohr_to_angstrom,
                float(cube[k][i+6].split()[3]) * bohr_to_angstrom,
                float(cube[k][i+6].split()[4]) * bohr_to_angstrom,
            ])
        structure.get_atoms(atoms_list)
        # end read structure info
        if k == 0:
            write_structure(structure=structure, filepath=args.output_structure+".total.cif")
        elif k == 1:
            write_structure(structure=structure, filepath=args.output_structure+".part1.cif")
        elif k == 2:
            write_structure(structure=structure, filepath=args.output_structure+".part2.cif")
    
    a = np.linalg.norm(structure.cell[0])
    b = np.linalg.norm(structure.cell[1])
    c = np.linalg.norm(structure.cell[2])
    
    # assume three cube file have the same ngridx ngridy and ngridz
    ngridx = int(cube[0][3].split()[0])
    ngridy = int(cube[0][4].split()[0])
    ngridz = int(cube[0][5].split()[0])          

    data_iii = []
    for i in range(3):
        # read grid value
        tmp_str = "".join(cube[i][natoms_each_image[i]+6:])
        data = np.fromstring(tmp_str, sep="\n")
            
        #data = data.reshape(ngridz, ngridy, ngridx)
        # grid data in cube is iterated in different compared to *CHG* of vasp
        data = data.reshape(ngridx, ngridy, ngridz)
        data_iii.append(data)
        # charge data in cube file is in shape (ngridx, ngridy, ngridz)
        # while charge in *CHG* file is in shape (ngzf, ngyf, ngxf)
        # they are different!        
    
    
    # data_iii dimension reduction
    data_sub = np.array(data_iii[0]) - np.array(data_iii[1]) - np.array(data_iii[2])
    # the unit of value is actually not physical now!
    #
    # cell_volume are in unit of Angstrom^3
    cell_volume = np.dot(np.cross(np.array(structure.cell[0]), np.array(structure.cell[1])), np.array(structure.cell[2]))
    cell_volume_per_unit = cell_volume / (ngridx * ngridy * ngridz)
    
    # value in cube file are \rho(r)_of_electrons in unit of e/Bohr^3
    # namely number of electrons each Borh^3
    # so we have to convert it to e/Angstrom^23, through divide it by borh_to_angstrom**3
    total_electrons = np.sum(data_iii[0])  * cell_volume_per_unit / bohr_to_angstrom**3
    #
    
    print("======================================================\n")
    print("           Information collected\n")
    print("------------------------------------------------------\n")
    print("cell volume: %f (A^3)\n" % cell_volume)
    print("total electrons: %f\n" % total_electrons)
    
    
    
    # data_sub is in unit of e/Bohr^3
    # we will build data_red_? to be in unit of e/Anstrom, namely number of electrons per Angstrom
    # to do this we have to time the volume density with bohr_to_angstrom^-3
    data_red_a = []
    data_red_b = []
    data_red_c = []
    if "c" in args.abscissa:
        len_ci = c / ngridz
        for ci in range(data_sub.shape[2]):
            tmp = 0
            for bi in range(data_sub.shape[1]):
                tmp += np.sum(data_sub[:, bi, ci])
            nelect_ci = tmp * cell_volume_per_unit / bohr_to_angstrom**3
            rho_line = nelect_ci / len_ci
            data_red_c.append(rho_line)
    if "b" in args.abscissa:
        len_bi = b / ngridy    
        for bi in range(data_sub.shape[1]):
            tmp = 0
            for ai in range(data_sub.shape[0]):
                tmp += np.sum(data_sub[ai, bi, :])
            nelect_bi = tmp * cell_volume_per_unit / bohr_to_angstrom**3
            rho_line = nelect_bi / len_bi                
            data_red_b.append(rho_line)
    if "a" in args.abscissa:
        len_ai = a / ngridx     
        for ai in range(data_sub.shape[0]):
            tmp = 0
            for ci in range(data_sub.shape[2]):
                tmp += np.sum(data_sub[ai, :, ci])
            nelect_ai = tmp * cell_volume_per_unit / bohr_to_angstrom**3
            rho_line = nelect_ai / len_ai                   
            data_red_a.append(rho_line)    

    # output the data and make the plot
    if "c" in args.abscissa:
        with open(args.output+".1d.c.data", 'w') as fout:
            fout.write("#c(angstrom) \Delta(rho) (number of electron per Angstrom)\n")
            c_coord = np.linspace(0, c, len(data_red_c))
            for i in range(len(data_red_c)):
                fout.write("%f %f\n" % (c_coord[i], data_red_c[i]))
        plt.plot(np.linspace(0, c, len(data_red_c)), data_red_c)                
        plt.ylabel(r"$\Delta\rho (e/\AA)$")
        plt.tight_layout()
        plt.savefig(args.output+".1d.c.png")
        plt.close()                
    if "b" in args.abscissa:
        with open(args.output+".1d.b.data", 'w') as fout:
            fout.write("#b(angstrom) \Delta(rho) (number of electron per Angstrom)\n")
            b_coord = np.linspace(0, b, len(data_red_b))
            for i in range(len(data_red_b)):
                fout.write("%f %f\n" % (b_coord[i], data_red_b[i]))        
        plt.plot(np.linspace(0, b, len(data_red_b)), data_red_b)    
        plt.ylabel(r"$\Delta\rho (e/\AA)$")              
        plt.tight_layout()
        plt.savefig(args.output+".1d.b.png")
        plt.close()                
    if "a" in args.abscissa:
        with open(args.output+".1d.a.data", 'w') as fout:
            fout.write("#a(angstrom) \Delta(rho) (number of electron per Angstrom)\n")
            a_coord = np.linspace(0, a, len(data_red_a))
            for i in range(len(data_red_a)):
                fout.write("%f %f\n" % (a_coord[i], data_red_a[i]))
        plt.plot(np.linspace(0, a, len(data_red_a)), data_red_a)                
        plt.ylabel(r"$\Delta\rho (e/\AA)$")            
        plt.tight_layout()
        plt.savefig(args.output+".1d.a.png")
        plt.close()
    


if __name__ == "__main__":
    main()