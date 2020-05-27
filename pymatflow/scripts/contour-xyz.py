#!/usr/bin/env python

import os
import sys
import argparse

from pymatflow.cmd.structflow import read_structure

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    parser.add_argument("-o", "--output", type=str, default="./contour.png",
        help="the output image file name")
    
    parser.add_argument("--atoms", nargs='+', type=int,
        help="list of fixed atoms, index start from 1, default is all atoms")

    parser.add_argument("--dgrid3d", type=int, nargs=3,
        default=[100, 100, 4],
        help="used by gnuplot to set dgrid3d int, int, int")

    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()


    structure = read_structure(args.input)

    if args.atoms == None:
        args.atoms = range(1, len(structure.atoms)+1)
    
    data = []
    for i in args.atoms:
        #X.append(structure.atoms[i-1].x)
        #Y.append(structure.atoms[i-1].y)
        #Z.append(structure.atoms[i-1].z)
        data.append([structure.atoms[i-1].x, structure.atoms[i-1].y, structure.atoms[i-1].z])

    with open(args.output+".data", "w") as fout:
        fout.write("# format: x y z\n")
        for d in data:
            fout.write("%f\t%f\t%f\n" % (d[0], d[1], d[2]))
    

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
        fout.write("set surface\n")
        fout.write("set pm3d hidden3d\n")
        fout.write("set output \"%s\"\n" % (args.output+".surf"+".gif"))
        fout.write("splot '%s' u 1:2:3 w pm3d\n" % (args.output+".data"))

        #fout.write("set origin 0.5, 0\n")
        fout.write("\n\n")
        fout.write("set pm3d map\n")
        fout.write("set surface\n")
        fout.write("set output \"%s\"\n" % (args.output+".map"+".gif"))
        fout.write("splot '%s' u 1:2:3 w pm3d\n" % (args.output+".data"))
    os.system("gnuplot plot.gnuplot")

if __name__ == "__main__":
    main()