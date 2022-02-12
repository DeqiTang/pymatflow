#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import re
import shutil
import argparse
from pymatflow.base.xyz import BaseXyz

"""
support for auto preparation of pseudopotential files for QE, Abinit, Siesta

support for qe not implemented now

args.qe_pot_dir:
    pot-qe/SSSP_efficiency_pseudos/*
    pot-qe/SSSP_precision_pseudos/*


args.abinit_pot_dir:
    pot-abinit/JTH-PBE-atomicdata-1.1/ATOMICDATA/
        elment.GGA_PBE-JTH.xml
    pot-abinit/pbe_s_sr/
        element.psp8

args.siesta_pot_dir:
    pot-siesta/abinit-ncpp/
        element.psf
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--xyz", type=str,
            help="the xyz file name")

    parser.add_argument("-p", "--program", type=str,
            choices=["qe", "abinit", "siesta"],
            help="choose to prepare pseudo for program")

    parser.add_argument("-d", "--directory", type=str, default="./",
            help="directory to put the pseudopotential files")

    parser.add_argument("--qe-type", type=str, default="paw_pbe",
            choices=["paw_pbe", "PAW_PBE", "sssp_efficiency", "SSSP_EFFICIENCY", "sssp_precision", "SSSP_PRECISION"],
            help="type of qe pseudopotential")

    parser.add_argument("--abinit-type", type=str, default="paw_pbe",
            choices=["paw_pbe", "ncpp", "PAW_PBE", "NCPP"],
            help="type of abinit pseudopotential")

    parser.add_argument("--siesta-type", type=str, default="abinit-ncpp",
            help="type of siesta pseudopotential")

    parser.add_argument("--qe-pot-dir", type=str,
            default=os.path.join(os.path.expanduser("~"), ".pot-qe"),
            help="vasp pot databse directory")

    parser.add_argument("--abinit-pot-dir", type=str,
            default=os.path.join(os.path.expanduser("~"), ".pot-abinit"),
            help="vasp pot databse directory")

    parser.add_argument("--siesta-pot-dir", type=str,
            default=os.path.join(os.path.expanduser("~"), ".pot-siesta"),
            help="vasp pot databse directory")


    args = parser.parse_args()

    print("=================================================================\n")
    print("                 pot-from-xyz-modified.py\n")
    print("-----------------------------------------------------------------\n")
    print("generate PseudoPotential file from information on the xyz file.\n")
    print("\n")

    xyz = BaseXyz()
    xyz.get_xyz(args.xyz)

    # --------------------------------------------------------------------------
    # deal with program --------> Quantum ESPRESSO
    # --------------------------------------------------------------------------
    if args.program.lower() == "qe":
        if args.qe_type.lower() == "paw_pbe":
            pass
        elif args.qe_type.lower() == "sssp_efficiency":
            # use SSSP
            # https://www.materialscloud.org/discover/sssp/table/efficiency
            all_pseudos = os.listdir(os.path.join(args.qe_pot_dir, "SSSP_efficiency_pseudos"))
            cmd = "cp "
            for element in xyz.specie_labels:
                for pseudo in all_pseudos:
                    #if re.match(element, pseudo, re.IGNORECASE):
                    #if re.match(element+".", pseudo, re.IGNORECASE) or re.match(element+"_", pseudo, re.IGNORECASE):
                    if pseudo.split(".")[0].lower() == element.lower() or pseudo.split("_")[0].lower() == element.lower():
                        cmd = cmd + "%s " % (os.path.join(args.qe_pot_dir, "SSSP_efficiency_pseudos/%s" % pseudo))
                        break
            cmd = cmd + " %s" % args.directory
            os.system(cmd)
        elif args.qe_type.lower() == "sssp_precision":
            # use SSSP
            # https://www.materialscloud.org/discover/sssp/table/precision
            all_pseudos = os.listdir(os.path.join(args.qe_pot_dir, "SSSP_precision_pseudos"))
            cmd = "cp "
            for element in xyz.specie_labels:
                for pseudo in all_pseudos:
                    #if re.match(element, pseudo, re.IGNORECASE):
                    #if re.match(element+".", pseudo, re.IGNORECASE) or re.match(element+"_", pseudo, re.IGNORECASE):
                    if pseudo.split(".")[0].lower() == element.lower() or pseudo.split("_")[0].lower() == element.lower():
                        cmd = cmd + "%s " % (os.path.join(args.qe_pot_dir, "SSSP_precision_pseudos/%s" % pseudo))
                        break
            cmd = cmd + " %s" % args.directory
            os.system(cmd)

    # --------------------------------------------------------------------------
    # deal with program --------> Abinit
    # --------------------------------------------------------------------------
    if args.program.lower() == "abinit":
        if args.abinit_type.lower() == "paw_pbe":
            cmd = "cp "
            for element in xyz.specie_labels:
                cmd = cmd + "%s/JTH-PBE-atomicdata-1.1/ATOMICDATA/%s.GGA_PBE-JTH.xml " % (args.abinit_pot_dir, element)
            cmd = cmd + " %s" % args.directory
            os.system(cmd)
        elif args.abinit_type.lower() == "ncpp":
            cmd = "cp "
            for element in xyz.specie_labels:
                cmd = cmd + "%s/pbe_s_sr/%s.psp8 " % (args.abinit_pot_dir, element)
            cmd = cmd + " %s" % args.directory
            os.system(cmd)
        else:
            pass
    # --------------------------------------------------------------------------
    # deal with program --------> Siesta
    # --------------------------------------------------------------------------
    if args.program.lower() == "siesta":
        if args.siesta_type.lower() == "abinit-ncpp":
            # use Abinit NCPP pseudopotential for SIESTA
            cmd = "cp "
            for element in xyz.specie_labels:
                cmd = cmd + "%s/abinit-ncpp/%s.psf " % (args.siesta_pot_dir, element)
            cmd = cmd + " %s" % args.directory
            os.system(cmd)
        else:
            pass

    # ==========================================================================
