#!/usr/bin/env python

import os
import sys
import copy
import argparse

import numpy as np
import matplotlib.pyplot as plt

from pymatflow.cmd.structflow import read_structure
from pymatflow.cmd.structflow import write_structure

from pymatflow.vasp.post.pdos import post_pdos

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str, required=True,
        help="input PARCHG file")

    parser.add_argument("--output-structure", type=str, default="parchg.cif",
        help="output stucture contained in PARCHG")

    parser.add_argument("-o", "--output", type=str, default="./contour-ldos",
        help="prefix of the output image file name")

    parser.add_argument("--dgrid3d", type=int, nargs=3,
        default=[100, 100, 4],
        help="used by gnuplot to set dgrid3d int, int, int")


    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()
    
    parchg_filepath = args.input
    
    with open(parchg_filepath, "r") as fin:
        parchg = fin.readlines()
        
    for i in range(len(parchg)):
        if len(parchg[i].split()) == 0:
            first_blank_line = i
            break
    os.system("mkdir -p /tmp/pymatflow/")
    with open("/tmp/pymatflow/POSCAR", "w") as fout:
        for i in range(first_blank_line):
            fout.write(parchg[i])
    structure = read_structure("/tmp/pymatflow/POSCAR")
    write_structure(structure=structure, filepath=args.output_structure)
    
    ngxf = int(parchg[first_blank_line+1].split()[0])
    ngyf = int(parchg[first_blank_line+1].split()[1])
    ngzf = int(parchg[first_blank_line+1].split()[2])    
    
    data = np.loadtxt(parchg[first_blank_line+2:])
    
    with open("parchg.data", "w") as fout:
        for i in range(first_blank_line+2, len(parchg)):
            fout.write(parchg[i])

    with open("plot.gnuplot", "w") as fout:
        fout.write("set term gif\n")
        fout.write("set dgrid3d %d, %d, %d\n" % (args.dgrid3d[0], args.dgrid3d[1], args.dgrid3d[2]))
        fout.write("set contour base\n")
        fout.write("unset key\n")
        fout.write("set autoscale\n")
        fout.write("set xlabel \"X\"\n")
        fout.write("set ylabel \"Y\"\n")
        fout.write("set zlabel \"Z\"\n")
        #fout.write("set output \"%s\"\n" % (args.output+".gif"))
        #fout.write("set multiplot layout 1, 2\n")

        #fout.write("set origin 0, 0\n")
        fout.write("\n\n")
        fout.write("\n\n")
        fout.write("set pm3d map\n")
        fout.write("set surface\n")
        fout.write("set output \"controur.gif\"\n")
        fout.write("splot \"parchg.data\" u 1:2:10 w pm3d")        
        #fout.write("splot '%s' u 1:2:3 w pm3d\n" % (args.output+".data"))        
    os.system("gnuplot plot.gnuplot")
    
if __name__ == "__main__":
    main()