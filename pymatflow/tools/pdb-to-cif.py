#!/usr/bin/env python

import os
import argparse


"""
rely on openbabel: trizen -S openbabel
openbabel command in archlinux: obabel
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--pdb", type=str,
            help="the pdb structure file")
    parser.add_argument("-o", "--cif", type=str,
            help="the output cif file")
    args = parser.parse_args()
 
    print("==================================\n")
    print("         pdb-to-cif.py\n")
    print("----------------------------------\n")
    print("using openbabel comamnd obabel to\n")
    print('convert from pdb to cif file\n')   

    os.system("obabel -ipdb %s -ocif -O%s" % (args.pdb, args.cif))
