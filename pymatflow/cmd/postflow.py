#!/usr/bin/env python

import os
import sys
import argparse



def get_kpath(kpath_manual=None, kpath_file=None):
    """
    :param kpath_manual: manual input kpath like --kpath '0.000000 0.000000 0.000000 GAMMA 5' '0.500000 0.000000 0.000000 X 5' '0.0000 0.000 0.50000 A |' '0.5 0.5 0.5 R '
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



def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="driver", title="subcommands", description="choose one and only one calculator")

    # --------------------------------------------------------------------------
    # Abinit
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("abinit", help="using abinit as calculator")

    subparser.add_argument("-r", "--runtype", type=int, default="static",
            choices=[0, 1, 2, 3, 4, 5, 6, 7],
            help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-opt; 3->hexagonal-opt; 4->tetragonal-opt; 5->dfpt-elastic-piezo-dielec; 6->dfpt-phonon; 7->phonopy")

    subparser.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath for band structure calculation")

    subparser.add_argument("--kpath-file", type=str,
            help="file to read the kpath for band structure calculation")


    subparser.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory to do the calculation")


    # --------------------------------------------------------------------------
    # CP2K
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("cp2k", help="using cp2k as calculator")

    subparser.add_argument("-r", "--runtype", type=int, default="static",
            choices=[0, 1, 2, 3, 4 ,5, 6],
            help="choices of runtype. 0->static_run; 1->geo-opt; 2->cell-opt; 3->cubic-cell; 4->hexagonal-cell; 5->tetragonal-cell; 6->phonopy")


    subparser.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory for the running.")


    subparser.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath for band structure calculation")

    subparser.add_argument("--kpath-file", type=str,
            help="file to read the kpath for band structure calculation")



    # --------------------------------------------------------------------------
    # Quantum ESPRESSO
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("qe", help="using quantum espresso as calculator")

    subparser.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1, 2, 3, 4, 5, 6, 7, 8],
            help="choices of runtype. 0->static_run; 1->relax; 2->vc-relax; 3->cubic-cell; 4->hexagonal-cell; 5->tetragonal-cell; 6->neb; 7->dfpt; 8->phonopy")


    subparser.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory for the running.")

    subparser.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath in crystal_b, like --kpath-manual '0.000000 0.000000 0.000000 GAMMA 5' '0.500000 0.000000 0.000000 X 5' '0.0000 0.000 0.50000 A |' '0.5 0.5 0.5 R '")

    subparser.add_argument("--kpath-file", type=str,
            help="manual input kpath in crystal_b read from the file")

    subparser.add_argument("--nebint", type=str, default="pwscf.int",
            help="xxx.int")

    subparser.add_argument("--nebdat", type=str, default="pwscf.dat",
            help="xxx.dat")

    subparser.add_argument("--inpname", type=str, default="min-energy-path.gp",
            help="inpname for the gnuplot script")

    subparser.add_argument("--md", type=str, default="neb-report.md",
            help="Markdown report file name")

    subparser.add_argument("--nebout", type=str, default="neb.out",
            help="output file of neb calculation")



    # --------------------------------------------------------------------------
    # SIESTA
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("siesta", help="using siesta as calculator")

    subparser.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1, 2, 3, 4, 5],
            help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-cell; 3->hexagonal-cell; 4->tetragonal-cell; 5->phonpy")

    subparser.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory for the running.")

    subparser.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath like --kpath-manual '0.000000 0.000000 0.000000 GAMMA 5' '0.500000 0.000000 0.000000 X 5' '0.0000 0.000 0.50000 A |' '0.5 0.5 0.5 R '")

    subparser.add_argument("--kpath-file", type=str,
            help="manual input kpath read from the file")


    # --------------------------------------------------------------------------
    # VASP
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("vasp", help="using vasp as calculator")

    subparser.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1, 2, 3, 4, 5],
            help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-cell; 3->hexagonal-cell; 4->tetragonal-cell;")

    subparser.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory for the running.")


    subparser.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath like --kpath-manual '0.000000 0.000000 0.000000 GAMMA 5' '0.500000 0.000000 0.000000 X 5' '0.0000 0.000 0.50000 A |' '0.5 0.5 0.5 R '")

    subparser.add_argument("--kpath-file", type=str,
            help="manual input kpath read from the file")




    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()

    # if no argument passed to matflow
    if len(sys.argv) == 1:
        # display help message when no args provided
        parser.print_help()
        sys.exit(1)



# ==============================================================================
# Abinit Abinit Abinit Abinit Abinit Abinit Abinit Abinit Abinit Abinit Abinit
# ==============================================================================
    if args.driver == "abinit":
        if args.runtype == 0:
            # static
            pass
        elif args.runtype == 1:
            # optimization
            from pymatflow.abinit.post.opt import opt
            post = opt()
            post.parse(os.path.join(args.directory, "optimization.out"))
            post.export(directory=args.directory)
        elif args.runtype == 5:
            # dfpt-elastic-piezo-dielec
            #from pymatflow.abinit.post.dfpt import dfpt_elastic_piezo_dielec_anaddb_out
            #post = dfpt_elastic_piezo_dielec_anaddb_out()
            os.system("post-abinit-dfpt-elastic-piezo-dielec.py -d %s" % args.directory)
        elif args.runtype == 6:
            # dfpt-phonon
            pass
        elif args.runtype == 7:
            # phonopy
            os.system("post-abinit-phonopy.py -d %s -f %s --qpath-file %s --supercell-n %d %d %d" % (args.directory, args.xyz, args.kpath_file, args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))
        else:
            pass

# ==============================================================================
# CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K
# ==============================================================================
    elif args.driver == "cp2k":
        if args.runtype == 0:
            pass
        elif args.runtype == 6 :
            # phonopy
            os.system("post-cp2k-phonopy.py -d %s -f %s --qpath-file %s --supercell-n %d %d %d" % (args.directory, args.xyz, args.kpath_file, args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))
        else:
            pass
# ==============================================================================
# Quantum ESPRESSO Qautnum ESPRESSO Quantum ESPRESSO Quantum ESPRESSO Quantum ESPRESSO
# ==============================================================================
    elif args.driver == "qe":
        if args.runtype == 0:
            from pymatflow.qe.post.scf import scf_out
            post = scf_out()

        elif args.runtype == 1:
            pass
        elif args.runtype == 2:
            pass
        elif args.runtype == 3:
            from pymatflow.qe.neb import neb_post
            os.chdir(args.directory)
            task = neb_post(nebout=args.nebout)
            os.chdir("../")
            task.export(directory=args.directory, nebint=args.nebint, nebdat=args.nebdat, md=args.md)

        elif args.runtype == 7:
            # dfpt phonon
            os.system("post-qe-matdyn.py -d %s --option gnuplot --freq 0 0.1" % args.directory)
        elif args.runtype == 8:
            # phonopy phonon
            os.system("post-qe-phonopy.py -d %s -f %s --qpath-file %s --supercell-n %d %d %d" % (args.directory, args.xyz, args.kpath_file, args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))
        else:
            pass
# ==============================================================================
# SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA
# ==============================================================================
    elif args.driver == "siesta":
        if args.runtype == 0:
            pass
        elif args.runtype == 5:
            # phonopy
            os.system("post-siesta-phonopy.py -d %s -f %s --qpath-file %s --supercell-n %d %d %d" % (args.directory, args.xyz, args.kpath_file, args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))

# ==============================================================================
# VASP
# ==============================================================================
    elif args.driver == "vasp":
        if args.runtype == 0:
            pass
        elif args.runtype == 2:
            # cubic cell
            os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))

    # --------------------------------------------------------------------------



if __name__ == "__main__":
    main()