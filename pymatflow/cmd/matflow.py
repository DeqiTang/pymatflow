#!/usr/bin/env python

import os
from pymatflow.cmd.matflow_dftbplus import dftbPlusDriver, dftbPlusSubparser
from pymatflow.cmd.matflow_vasp import vaspDriver, vaspSubparser
from pymatflow.cmd.matflow_siesta import siestaDriver, siestaSubparser
from pymatflow.cmd.matflow_qe import qeDriver, qeSubparser
from pymatflow.cmd.matflow_cp2k import cp2kDriver, cp2kSubparser
from pymatflow.cmd.matflow_abinit import abinitDriver, abinitSubparser
import sys
import argparse


def get_kpath(kpath_manual=None, kpath_file=None):
    """
    :param kpath_manual: manual input kpath like --kpath-manual '0.000000 0.000000 0.000000 GAMMA 5' '0.500000 0.000000 0.000000 X 5' '0.0000 0.000 0.50000 A |' '0.5 0.5 0.5 R '
    :param kpath_file: manual input kpath read from the file
    :return kpath or None(when kpath_manual and kpath_file are both None)
    """
    # dealing with standard kpath
    kpath = None
    if kpath_manual != None:
        # kpath from script argument args.kpath
        kpath = []
        for kpoint in kpath_manual:
            if kpoint.split()[4] != "|":
                kpath.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    int(kpoint.split()[4]),
                    ])
            elif kpoint.split()[4] == "|":
                kpath.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    "|",
                    ])
        return kpath
    elif kpath_file != None:
        # kpath read from file specified by kpath_file
        # file is in format like this
        """
        5
        0.0 0.0 0.0 #GAMMA 15
        x.x x.x x.x #XXX |
        x.x x.x x.x #XXX 10
        x.x x.x x.x #XXX 15
        x.x x.x x.x #XXX 20
        """
        # if there is a '|' behind the label it means the path is
        # broken after that point!!!
        kpath = []
        with open(kpath_file, 'r') as fin:
            lines = fin.readlines()
        nk = int(lines[0])
        for i in range(nk):
            if lines[i+1].split("\n")[0].split()[4] != "|":
                kpath.append([
                    float(lines[i+1].split()[0]),
                    float(lines[i+1].split()[1]),
                    float(lines[i+1].split()[2]),
                    lines[i+1].split()[3].split("#")[1].upper(),
                    int(lines[i+1].split()[4]),
                    ])
            elif lines[i+1].split("\n")[0].split()[4] == "|":
                kpath.append([
                    float(lines[i+1].split()[0]),
                    float(lines[i+1].split()[1]),
                    float(lines[i+1].split()[2]),
                    lines[i+1].split()[3].split("#")[1].upper(),
                    '|',
                    ])
        return kpath
    else:
        pass
    # -------------------------------------------------------------------


def getXyzFile(args):
    """
    return-> xyzfile, images
    """
    # dealing wich structure files
    xyzfile = None
    images = []
    if args.xyz != None:
        xyzfile = args.xyz
    elif args.cif != None:
        os.system("structflow convert -i %s -o %s.xyz" % (args.cif, args.cif))
        xyzfile = "%s.xyz" % args.cif
    elif args.xsd != None:
        os.system("structflow convert -i %s -o %s.xyz" % (args.xsd, args.xsd))
        xyzfile = "%s.xyz" % args.xsd
    elif args.xsf != None:
        os.sytem("structflow convert -i % -o %s.xyz" % (args.xsf, args.xsf))
        xyzfile = "%s.xyz" % args.xsf
    else:
        # neb caculattion with
        images = []
        for image in args.images:
            if image.split(".")[-1] == "xyz":
                images.append(image)
            else:
                os.system("structflow convert -i %s -o %s.xyz" % (image, image))
                images.append("%s.xyz" % image)
        xyzfile = images[0] # this set only for dealing with pseudo potential file
        #
    # dealing with pseudo potential file
    if args.driver in ["dftb+"]:
        pass
    else:
        if args.pot == "./":
            #TODO make a simple check, whether there exists the potential file
            pass
        elif args.pot == "auto":
            if args.driver == "abinit":
                os.system("pot-from-xyz-modified.py -i %s -d ./ -p abinit --abinit-type=ncpp" % xyzfile)
            elif args.driver == "qe":
                if args.runtype == 6:
                    os.system("pot-from-xyz-modified.py -i %s -d ./ -p qe --qe-type=%s" % (images[0]), args.pot_type)
                else:
                    os.system("pot-from-xyz-modified.py -i %s -d ./ -p qe --qe-type=%s" % (xyzfile, args.pot_type))
            elif args.driver == "siesta":
                print("=============================================================\n")
                print("                     WARNING\n")
                print("-------------------------------------------------------------\n")
                print("support for auto preparation of pseudopotential file for siesta\n")
                print("is not fully implemented now!\n")
                print("please prepare it yourself\n")
                sys.exit(1)
            elif args.driver == "vasp":
                os.system("vasp-potcar-from-xyz.py --type %s -i %s -o ./POTCAR" % (args.pot_type, xyzfile))
        else:
            os.system("cp %s/* ./" % args.pot)
    return xyzfile, images


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="driver", title="subcommands", description="choose one and only one calculator")

    abinitSubparser(subparsers)
    cp2kSubparser(subparsers)
    qeSubparser(subparsers)
    siestaSubparser(subparsers)
    vaspSubparser(subparsers)
    dftbPlusSubparser(subparsers)

    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()

    # if no argument passed to matflow
    if len(sys.argv) == 1:
        # display help message when no args provided
        parser.print_help()
        sys.exit(1)

    # check directory
    if os.getcwd() == os.path.abspath(args.directory):
        print("================================================================\n")
        print("                    Warning from matflow\n")
        print("----------------------------------------------------------------\n")
        print("you are trying to specify current directory for task running\n")
        print("which is not allowed.\n")
        print("try specify it to a new directory like -d matflow-running\n")
        sys.exit(1)
        

    if args.driver == "abinit":
        abinitDriver(args)
    elif args.driver == "cp2k":
        cp2kDriver(args)
    elif args.driver == "qe":
        qeDriver(args)
    elif args.driver == "siesta":
        siestaDriver(args)
    elif args.driver == "vasp":
        vaspDriver(args)
    elif args.driver == "dftb+":
        dftbPlusDriver(args)
    # --------------------------------------------------------------------------


if __name__ == "__main__":
    main()
