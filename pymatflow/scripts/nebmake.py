#!/usr/bin/env python

import os
import argparse

from pymatflow.structure.neb import interpolate
from pymatflow.cmd.structflow import read_structure
from pymatflow.cmd.structflow import write_structure



def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--images", type=str, nargs=2,
            required=True,
            help="the initial and final structure file")

    parser.add_argument("-n", "--nimage", type=int, default=None,
            required=True,
            help="number of inter images")

    parser.add_argument("-m", "--moving-atom", type=int, nargs="+",
            required=True,
            help="specifying the moving atoms, index start from 0")

    parser.add_argument("-d", "--directory", type=str, default="./",
            help="directory to put the generated images")
    
    parser.add_argument("--frac", type=int, default=1,
            choices=[0, 1],
            help="1(default): use faractional, 0: use cartesian")
            
    # ==============================================================
    args = parser.parse_args()

    initial = read_structure(args.images[0])
    final = read_structure(args.images[1])

    inter_images = interpolate(initial=initial, final=final, nimage=args.nimage, moving_atom=args.moving_atom)

    
    os.system("mkdir -p %s" % os.path.join(args.directory, "%.2d" % (0)))
    write_structure(structure=initial, filepath=os.path.join(args.directory, "%.2d/POSCAR" % (0)), frac=args.frac)
    os.system("mkdir -p %s" % os.path.join(args.directory, "%.2d" % (args.nimage+1)))    
    write_structure(structure=final, filepath=os.path.join(args.directory, "%.2d/POSCAR" % (args.nimage+1)), frac=args.frac)
    for i in range(len(inter_images)):
        os.system("mkdir -p %s" % os.path.join(args.directory, "%.2d" % (i+1)))        
        write_structure(structure=inter_images[i], filepath=os.path.join(args.directory, "%.2d/POSCAR" % (i+1)), frac=args.frac)
    
    print("===========================================\n")
    print("generate inter images for neb calculation\n")
    print("===========================================\n")

    print("-------------------------------------------\n")

if __name__ == "__main__":
    main()