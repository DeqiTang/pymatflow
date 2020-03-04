#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import argparse
from vaspstudio.base.poscar import vasp_poscar

"""
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--xyz", type=str,
            help="the xyz file name")

    parser.add_argument("-o", "--poscar", type=str, default="POSCAR",
            help="the poscar file name, default is POSCAR")

    args = parser.parse_args()

    print("========================================\n")
    print("           xyz-modified-to-poscar.py\n")
    print("----------------------------------------\n")
    print("convert from xyz-modifed to poscar\n")
    print("with the help from vaspstudio\n")

    poscar = vasp_poscar()
    poscar.xyz.get_xyz(args.xyz)

    with open(args.poscar, 'w') as fout:
        poscar.to_poscar(fout)
