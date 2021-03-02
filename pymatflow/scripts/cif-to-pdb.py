#!/usr/bin/env python

import os
import argparse


"""
rely on openbabel: trizen -S openbabel
openbabel command in archlinux: obabel
Reference: http://openbabel.org/wiki/Babel#File_Formats
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--cif", type=str,
            help="the cif structure file")

    parser.add_argument("-o", "--pdb", type=str,
            help="the output pdb file")

    args = parser.parse_args()

    print("========================================\n")
    print("         cif-to-pdb.py\n")
    print("----------------------------------------\n")
    print("using openbabel comamnd obabel to\n")
    print('convert from cif to pdb file\n')
    print("Note: vest cannot visualize the output\n")
    print("pdb file normally, you can use other\n")
    print("softwares like xcrysden to visualize it\n")

    os.system("obabel -icif %s -opdb -O%s" % (args.cif, args.pdb))
