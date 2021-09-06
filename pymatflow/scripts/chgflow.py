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

from pymatflow.charge.chg_vasp import VaspCHG



def main():
    parser = argparse.ArgumentParser()
    
    subparsers = parser.add_subparsers(dest="driver", title="subcommands", description="choose one and only one subcommand")
    
    subparser = subparsers.add_parser("1d", help="dimension reduction of CHG* to one dimensional")

    subparser.add_argument("-i", "--input", type=str, required=True,
        help="input vasp *CHG* file, -i *CHG* ")

    subparser.add_argument("--output-structure", type=str, default="chg",
        help="output stucture contained in *CHG*")

    subparser.add_argument("-o", "--output", type=str, default="chg",
        help="prefix of the output image file name")
    
    subparser.add_argument("--levels", type=int, default=10,
        help="levels of the color map or color bar")

    subparser.add_argument("-z", "--z", type=float, default=1,
        help="a value between 0 and 1, indicat height of in z direction to print the plot")
        
    subparser.add_argument("--cmap", type=str, default="gray",
        choices=["gray", "hot", "afmhot", "Spectral", "plasma", "magma", "hsv", "rainbow", "brg"])
        
    subparser.add_argument("--abscissa", type=str, nargs="+", default=["a", "b", "c"], 
        choices=["a", "b", "c"], 
        help="choose the direction to do the dimension reduction")

    #  add 
    subparser = subparsers.add_parser("add", help="add CHG*")

    subparser.add_argument("-i", "--input", type=str, nargs="+", required=True,
        help="input vasp *CHG* file, -i CHG1 CHG2 CHG3 ... ")

    subparser.add_argument("--output-structure", type=str, default="chg",
        help="output stucture contained in *CHG*")

    subparser.add_argument("-o", "--output", type=str, default="chg-merged",
        help="prefix of the output chg file name")

    # slice
    subparser = subparsers.add_parser("slice", help="slice CHG* to get 2d DATA")

    subparser.add_argument("-i", "--input", type=str, required=True,
        help="input vasp *CHG* file, -i *CHG* ")

    subparser.add_argument("--output-structure", type=str, default="chg-slice",
        help="output stucture contained in *CHG*")

    subparser.add_argument("-o", "--output", type=str, default="chg",
        help="prefix of the output image file name")
    
    subparser.add_argument("--levels", type=int, default=10,
        help="levels of the color map or color bar")

    subparser.add_argument("-z", "--z", type=float, default=1,
        help="a value between 0 and 1, indicat height of in z direction to print the plot")
        
    subparser.add_argument("--cmap", type=str, default="gray",
        choices=["gray", "hot", "afmhot", "Spectral", "plasma", "magma", "hsv", "rainbow", "brg"])
        
    subparser.add_argument("--abscissa", type=str, nargs="+", default=["a", "b", "c"], 
        choices=["a", "b", "c"], 
        help="choose the direction to do the dimension reduction")
    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()
    
    if args.driver == "1d":

        vaspchg = VaspCHG(args.input)
                    

        a = np.linalg.norm(vaspchg.structure.cell[0])
        b = np.linalg.norm(vaspchg.structure.cell[1])
        c = np.linalg.norm(vaspchg.structure.cell[2])
        
        
        print("======================================================\n")
        print("           Information collected\n")
        print("------------------------------------------------------\n")
        print("cell volume: %f (A^3)\n" % vaspchg.cell_volume)
        print("total electrons: %f\n" % vaspchg.total_electrons)
        
        
        # unit of data_red_? is e/Anstrom, namely number of electrons per Angstrom
        data_red_a = []
        data_red_b = []
        data_red_c = []
        if "c" in args.abscissa:
            factor = vaspchg.cell_volume_per_unit / vaspchg.cell_volume
            len_ci = c / vaspchg.ngzf
            for ci in range(vaspchg.data.shape[0]):
                tmp = 0
                for bi in range(vaspchg.data.shape[1]):
                    tmp += np.sum(vaspchg.data[ci, bi, :])
                nelect_ci = tmp * factor
                rho_line = nelect_ci / len_ci
                data_red_c.append(rho_line)
        if "b" in args.abscissa:
            factor = vaspchg.cell_volume_per_unit / vapchg.cell_volume
            len_bi = b / vaspchg.ngyf    
            for bi in range(vaspchg.data.shape[1]):
                tmp = 0
                for ai in range(vaspchg.data.shape[2]):
                    tmp += np.sum(vaspchg.data[:, bi, ai])
                nelect_bi = tmp * factor
                rho_line = nelect_bi / len_bi                
                data_red_b.append(rho_line)
        if "a" in args.abscissa:
            factor = vaspchg.cell_volume_per_unit / vaspchg.cell_volume
            len_ai = a / vaspchg.ngxf        
            for ai in range(vaspchg.data.shape[2]):
                tmp = 0
                for ci in range(vaspchg.data.shape[0]):
                    tmp += np.sum(vaspchg.data[ci, :, ai])
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

    elif args.driver == "add":
        chgs = []   
        for item in args.input:
            chgs.append(VaspCHG(item))
        pass
       
    elif args.driver == "slice":
        vaspchg = VaspCHG(args.input)
        vaspchg.plot_grayscale_z(z=args.z, output_prefix=args.output)
        vaspchg.plot_contour_2d(z=args.z, levels=args.levels, cmap=args.cmap, output_prefx=args.output)


if __name__ == "__main__":
    main()
