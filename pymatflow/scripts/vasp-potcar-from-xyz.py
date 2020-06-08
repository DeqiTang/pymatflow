#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import argparse
from pymatflow.vasp.base.poscar import vasp_poscar

"""
currently only support auto generation of PAW_PBE type POTCAR
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--xyz", type=str,
            help="the xyz file name")

    parser.add_argument("-o", "--potcar", type=str,
            default="POTCAR",
            help="the cif file name")

    parser.add_argument("--pot-database", type=str,
            default=os.path.join(os.path.expanduser("~"), ".pot-vasp"),
            help="vasp pot databse directory")

    parser.add_argument("--type", type=str, default="PAW_PBE",
            choices=["PAW_PBE", "PAW_LDA", "PAW_PW91", "paw_pbe", "paw_lda", "paw_pw91"],
            help="choose type of POT for POTCAR")

    args = parser.parse_args()

    print("================================================\n")
    print("             potcar-from-xyz.py\n")
    print("------------------------------------------------\n")
    print("generate POTCAR from information on the xyz file\n")
    print("\n")
    print("Note: elements are in the order of atomic number\n")

    poscar = vasp_poscar()
    poscar.xyz.get_xyz(args.xyz)

    cmd = "cat "
    for element in poscar.xyz.specie_labels:
        if args.type.upper() == "PAW_PBE":
            cmd = cmd + os.path.join(args.pot_database, "PAW_PBE/%s/POTCAR " % element)
        elif args.type.upper() == "PAW_LDA":
            cmd = cmd + os.path.join(args.pot_database, "PAW_LDA/%s/POTCAR " % element)
        elif args.type.upper() == "PAW_PW91":
            cmd = cmd + os.path.join(args.pot_database, "PAW_PW91/%s/POTCAR " % element)
        else:
            pass
    cmd = cmd + "> %s" % args.potcar
    os.system(cmd)
