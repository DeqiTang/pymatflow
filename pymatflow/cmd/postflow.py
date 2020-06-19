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

    subparser.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1, 2, 3, 4, 5, 6, 7],
            help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-opt; 3->hexagonal-opt; 4->tetragonal-opt; 5->dfpt-elastic-piezo-dielec; 6->dfpt-phonon; 7->phonopy")

    subparser.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath for band structure calculation")

    subparser.add_argument("--kpath-file", type=str,
            help="file to read the kpath for band structure calculation")


    subparser.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory to do the calculation")

    # structure file
    structfile = subparser.add_mutually_exclusive_group() # only one of them can be provided
    structfile.add_argument("--xyz", type=str, default=None,
            help="The xyz structure file with the second line specifying the cell parameter")

    structfile.add_argument("--cif", type=str, default=None,
            help="The cif structure file")

    structfile.add_argument("--xsd", type=str, default=None,
            help="The xsd structure file")

    structfile.add_argument("--xsf", type=str, default=None,
            help="The xsf structure file")


    subparser.add_argument("--plotrange", type=float, nargs="+",
            default=[0, 1.0],
            help="range to plot. in percentage")

    subparser.add_argument("--bandrange", type=float, nargs="+",
            default=[0, 1.0],
            help="band range to plot. in percentage")

    subparser.add_argument("--engine", type=str, default="matplotlib",
            choices=["matplotlib", "gnuplot"],
            help="plot engine, matplotlib or gnuplot")

    subparser.add_argument("--supercell-n", type=int, nargs="+",
            default=[1, 1, 1],
            help="supercell for phonopy, like [2, 2, 2]")

    # --------------------------------------------------------------------------
    # CP2K
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("cp2k", help="using cp2k as calculator")

    subparser.add_argument("-r", "--runtype", type=int, default="static",
            choices=[0, 1, 2, 3, 4 ,5, 6, 7, 8, 9, 10],
            help="choices of runtype. 0->static_run; 1->geo-opt; 2->cell-opt; 3->cubic-cell; 4->hexagonal-cell; 5->tetragonal-cell; 6->neb; 7->phonopy; 8->vibrational_analysis; 9:converge test; 10->aimd")

    subparser.add_argument("--static", type=str, default="scf",
            choices=["scf", "band"],
            help="type of static calc, like band")

    subparser.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory for the running.")


    subparser.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath for band structure calculation")

    subparser.add_argument("--kpath-file", type=str,
            help="file to read the kpath for band structure calculation")

    subparser.add_argument("--bandrange", type=float, nargs="+",
            default=[0, 1.0],
            help="band range to plot. in percentage")

    subparser.add_argument("--xrange", type=float, nargs=2,
            default=None,
            help="x range of the plot")

    subparser.add_argument("--yrange", type=float, nargs=2,
            default=None,
            help="y range of the plot")

    subparser.add_argument("--engine", type=str, default="matplotlib",
            choices=["matplotlib", "gnuplot"],
            help="plot engine, matplotlib or gnuplot")

    # structure file
    structfile = subparser.add_mutually_exclusive_group() # only one of them can be provided
    structfile.add_argument("--xyz", type=str, default=None,
            help="The xyz structure file with the second line specifying the cell parameter")

    structfile.add_argument("--cif", type=str, default=None,
            help="The cif structure file")

    structfile.add_argument("--xsd", type=str, default=None,
            help="The xsd structure file")

    structfile.add_argument("--xsf", type=str, default=None,
            help="The xsf structure file")

    subparser.add_argument("--converge", type=str, default="cutoff",
            choices=["cutoff", "rel_cutoff", "kpoints_auto", "kpoints_manual"],
            help="choose type of converge test")

    subparser.add_argument("--criteria", type=float, default=7.35e-4,
            help="converge criteria for cutoff or rel_cutoff or kpoints in unit of Ry")

    subparser.add_argument("--supercell-n", type=int, nargs="+",
            default=[1, 1, 1],
            help="supercell for phonopy, like [2, 2, 2]")


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

    subparser.add_argument("--bandrange", type=float, nargs="+",
            default=[0, 1.0],
            help="band range to plot. in percentage")

    subparser.add_argument("--engine", type=str, default="matplotlib",
            choices=["matplotlib", "gnuplot"],
            help="plot engine, matplotlib or gnuplot")


    subparser.add_argument("--plotrange", type=float, nargs="+",
            default=[0, 1.0],
            help="plot range (in percentage), like --plotrange 0.1 0.9")

    subparser.add_argument("--atomtoproj", type=int, nargs="+",
            default=[],
            help="atom to projection in atom projected dos. atom number starts with 1.")

    subparser.add_argument("--fontsize", type=int, default=10,
            help="fontsize for the plot.")
            
    # structure file
    structfile = subparser.add_mutually_exclusive_group() # only one of them can be provided
    structfile.add_argument("--xyz", type=str, default=None,
            help="The xyz structure file with the second line specifying the cell parameter")

    structfile.add_argument("--cif", type=str, default=None,
            help="The cif structure file")

    structfile.add_argument("--xsd", type=str, default=None,
            help="The xsd structure file")

    structfile.add_argument("--xsf", type=str, default=None,
            help="The xsf structure file")


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

    subparser.add_argument("--supercell-n", type=int, nargs="+",
            default=[1, 1, 1],
            help="supercell for phonopy, like [2, 2, 2]")

    subparser.add_argument("--opt-out", type=str, default=None,
            help="path of the output file for relax or vc-relax")


    subparser.add_argument("--xrange", type=float, nargs=2,
            default=None,
            help="x range of the plot")

    subparser.add_argument("--yrange", type=float, nargs=2,
            default=None,
            help="y range of the plot")
            
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


    # structure file
    structfile = subparser.add_mutually_exclusive_group() # only one of them can be provided
    structfile.add_argument("--xyz", type=str, default=None,
            help="The xyz structure file with the second line specifying the cell parameter")

    structfile.add_argument("--cif", type=str, default=None,
            help="The cif structure file")

    structfile.add_argument("--xsd", type=str, default=None,
            help="The xsd structure file")

    structfile.add_argument("--xsf", type=str, default=None,
            help="The xsf structure file")

    subparser.add_argument("--supercell-n", type=int, nargs="+",
            default=[1, 1, 1],
            help="supercell for phonopy, like [2, 2, 2]")

    subparser.add_argument("--engine", type=str, default="gnuplot",
            choices=["matplotlib", "gnuplot"],
            help="choose gnuplot or matplot lib to do the band plot")
    # --------------------------------------------------------------------------
    # VASP
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("vasp", help="using vasp as calculator")

    subparser.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1, 2, 3, 4, 5, 6, 7],
            help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-cell; 3->hexagonal-cell; 4->tetragonal-cell; 5->neb; 6->vasp-phonon; 7->phonopy")

    subparser.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory for the running.")

    
    subparser.add_argument("--static", type=str, default="band",
            choices=["scf", "band", "dos", "optics", "bse"],
            help="in case of band(default), run scf, nscf(bands) in a single run; in case of scf, run scf only, in case of optics, run scf and optics calc in a single run")


    subparser.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath like --kpath-manual '0.000000 0.000000 0.000000 GAMMA 5' '0.500000 0.000000 0.000000 X 5' '0.0000 0.000 0.50000 A |' '0.5 0.5 0.5 R '")

    subparser.add_argument("--kpath-file", type=str,
            help="manual input kpath read from the file")

    # structure file
    structfile = subparser.add_mutually_exclusive_group() # only one of them can be provided
    structfile.add_argument("--xyz", type=str, default=None,
            help="The xyz structure file with the second line specifying the cell parameter")

    structfile.add_argument("--cif", type=str, default=None,
            help="The cif structure file")

    structfile.add_argument("--xsd", type=str, default=None,
            help="The xsd structure file")

    structfile.add_argument("--xsf", type=str, default=None,
            help="The xsf structure file")


    subparser.add_argument("--efermi", type=str, default="nscf",
            choices=["scf", "nscf"],
            help="choose to read the efermi from nscf or scf")

    subparser.add_argument("--plotrange", type=float, nargs="+",
            default=[0, 1.0],
            help="range to plot. in percentage")

    subparser.add_argument("--bandrange", type=float, nargs="+",
            default=[0, 1.0],
            help="band range to plot. in percentage")

    subparser.add_argument("--xrange", type=float, nargs=2,
            default=None,
            help="x range of the plot")

    subparser.add_argument("--yrange", type=float, nargs=2,
            default=None,
            help="y range of the plot")

    subparser.add_argument("--engine", type=str, default="matplotlib",
            choices=["matplotlib", "gnuplot"],
            help="plot engine, matplotlib or gnuplot")

    subparser.add_argument("--supercell-n", type=int, nargs="+",
            default=[1, 1, 1],
            help="supercell for phonopy, like [2, 2, 2]")
            

    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()

    # if no argument passed to matflow
    if len(sys.argv) == 1:
        # display help message when no args provided
        parser.print_help()
        sys.exit(1)

    # dealing wich structure files
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




# ==============================================================================
# Abinit Abinit Abinit Abinit Abinit Abinit Abinit Abinit Abinit Abinit Abinit
# ==============================================================================
    if args.driver == "abinit":
        if args.runtype == 0:
            # static
            from pymatflow.abinit.post.bands import post_bands
            from pymatflow.base.xyz import base_xyz
            xyz = base_xyz()
            xyz.get_xyz(xyzfile)
            post = post_bands()
            post.get_xcoord_k(kpath=get_kpath(args.kpath_manual, args.kpath_file), cell=xyz.cell)
            post.get_ebands_agr(filepath=os.path.join(args.directory, "static-o_DS3_EBANDS.agr"))
            post.export(directory=args.directory, bandrange=args.bandrange, option=args.engine)
        elif args.runtype == 1:
            # optimization
            from pymatflow.abinit.post.opt import opt
            post = opt()
            post.parse(os.path.join(args.directory, "optimization.out"))
            post.export(directory=args.directory)
        elif args.runtype == 2:
            # cubic cell
            #os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
            os.system("cd %s; bash get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
        elif args.runtype == 3:
            # hexagonal cell
            #os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
            os.system("cd %s; bash get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
        elif args.runtype == 4:
            # tetragonal cell
            #os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
            os.system("cd %s; bash get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
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
            os.system("post-abinit-phonopy.py -d %s -f %s --qpath-file %s --supercell-n %d %d %d" % (args.directory, xyzfile, args.kpath_file, args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))
        else:
            pass

# ==============================================================================
# CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K
# ==============================================================================
    elif args.driver == "cp2k":
        if args.runtype == 0:
            # static
            if args.static == "scf":
                from pymatflow.cp2k.post.scf import scf_out
                task = scf_out()
                task.get_info(os.path.join(args.directory, "static-scf.out"))
                task.export(args.directory)
            elif args.static == "band":
                from pymatflow.cp2k.post.bands import bands_post
                from pymatflow.cmd.structflow import read_structure
                structure = read_structure(xyzfile)
                task = bands_post()
                task.get_efermi(static_out=os.path.join(args.directory, "static-scf.out"))
                task.get_kpath_and_bands(kpath=get_kpath(kpath_manual=args.kpath_manual, kpath_file=args.kpath_file), cell=structure.cell, bands=os.path.join(args.directory, "bands.bs"))
                task.export(directory=args.directory, engine=args.engine, bandrange=args.bandrange, xrange=args.xrange, yrange=args.yrange)
                task.print_gap()
                task.print_effective_mass()
        elif args.runtype == 1:
            from pymatflow.cp2k.post.opt import opt_out 
            task = opt_out()
            task.get_info(os.path.join(args.directory, "geo-opt.out"))
            task.export(args.directory)
        elif args.runtype == 2:
            from pymatflow.cp2k.post.opt import opt_out 
            task = opt_out()
            task.get_info(os.path.join(args.directory, "cell-opt.out"))
            task.export(args.directory)
        elif args.runtype == 3:
            # cubic cell
            os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
        elif args.runtype == 4:
            # hexagonal cell
            os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
        elif args.runtype == 5:
            # tetragonal cell
            os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
        elif args.runtype == 6:
            # neb 
            pass
        elif args.runtype == 7:
            # phonopy
            os.system("post-cp2k-phonopy.py -d %s -f %s --qpath-file %s --supercell-n %d %d %d" % (args.directory, xyzfile, args.kpath_file, args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))
        elif args.runtype == 8:
            # vibrational analysis
            pass
        elif args.runtype == 9:
            # converge test
            from pymatflow.cp2k.post.converge import converge_post
            task = converge_post()
            task.criteria_for_cutoff = args.criteria
            task.criteria_for_rel_cutoff = args.criteria
            task.criteria_for_kpoints = args.criteria
            task.postprocess(directory=args.directory, converge=args.converge)
        elif args.runtype == 10:
            # aimd
            from pymatflow.cp2k.post.md import md_post
            task = md_post(output=os.path.join(args.directory, "aimd.out"), run_type='MD')
            task.export(args.directory)
        else:
            pass
# ====================================================================================
# Quantum ESPRESSO Qautnum ESPRESSO Quantum ESPRESSO Quantum ESPRESSO Quantum ESPRESSO
# ====================================================================================
    elif args.driver == "qe":
        if args.runtype == 0:
            from pymatflow.qe.post.scf import scf_out
            from pymatflow.qe.post.bands import bands_post
            from pymatflow.qe.post.pdos import pdos_out
            os.chdir(args.directory)
            task = bands_post(pwxbandsin="static-bands.in", bandsxout="bands.out")
            task.plot_band(option=args.engine, bandrange=args.bandrange, xrange=args.xrange, yrange=args.yrange)
            os.chdir("../")
            task = pdos_out()
            task.get_data(directory=args.directory, filpdos="projwfc")
            task.export(directory=args.directory, plotrange=args.plotrange, atomtoproj=args.atomtoproj, fontsize=args.fontsize)

        elif args.runtype == 1:
            # relax
            from pymatflow.qe.post.opt import opt_out
            task = opt_out()
            #task.get_info(os.path.join(args.directory, "relax.out"))
            if args.opt_out == None:
                # use default relax.out
                args.opt_out = os.path.join(args.directory, "relax.out")
            task.get_info(args.opt_out)
            task.export(args.directory)
        elif args.runtype == 2:
            # vc-relax
            from pymatflow.qe.post.opt import opt_out
            task = opt_out()
            #task.get_info(os.path.join(args.directory, "vc-relax.out"))
            if args.opt_out == None:
                # use default vc-relax.out
                args.opt_out = os.path.join(args.directory, "vc-relax.out")
            task.get_info(args.opt_out)
            task.export(args.directory)
        elif args.runtype == 3:
            # cubic cell
            os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))    
        elif args.runtype == 4:
            # hexagonal cell
            os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
        elif args.runtype == 5:
            # tetragonal cell
            os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
        elif args.runtype == 6:
            # nudged elastic band
            from pymatflow.qe.post.neb import neb_post
            os.chdir(args.directory)
            task = neb_post(nebout=args.nebout)
            os.chdir("../")
            task.export(directory=args.directory, nebint=args.nebint, nebdat=args.nebdat, md=args.md)
        elif args.runtype == 7:
            # dfpt phonon
            os.system("post-qe-matdyn.py -d %s --option gnuplot --freq 0 0.1" % args.directory)
        elif args.runtype == 8:
            # phonopy phonon
            os.system("post-qe-phonopy.py -d %s -f %s --qpath-file %s --supercell-n %d %d %d" % (args.directory, xyzfile, args.kpath_file, args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))
        else:
            pass
# ==============================================================================
# SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA
# ==============================================================================
    elif args.driver == "siesta":
        if args.runtype == 0:
            # static
            from pymatflow.siesta.post.pdos import pdos
            task = pdos()
            task.get_info(os.path.join(args.directory, "siesta.PDOS.xml"))
            task.export(directory=args.directory)
            from pymatflow.siesta.post.bands import bands_post
            task = bands_post()
            task.process(os.path.join(args.directory, "siesta.bands"))
            task.export(directory=args.directory, option=args.engine)
        elif args.runtype == 1:
            # optimization
            from pymatflow.siesta.post.opt import opt_out
            task = opt_out()
            task.get_info(os.path.join(args.directory, "geometric-optimization.out"))
            task.export(directory=args.directory)
        elif args.runtype == 2:
            # cubic cell
            os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
        elif args.runtype == 3:
            # hexagonal cell
            os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
        elif args.runtype == 4:
            # tetragonal
            os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
        elif args.runtype == 5:
            # phonopy
            os.system("post-siesta-phonopy.py -d %s -f %s --qpath-file %s --supercell-n %d %d %d" % (args.directory, xyzfile, args.kpath_file, args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))
        else:
            pass

# ==============================================================================
# VASP
# ==============================================================================
    elif args.driver == "vasp":
        if args.runtype == 0:
            # static
            if args.static == "band":
                from pymatflow.vasp.post.bands import post_bands
                bands = post_bands()
                bands.get_kpath_and_vasprun(kpath=get_kpath(kpath_manual=args.kpath_manual, kpath_file=args.kpath_file), vasprun=os.path.join(args.directory, "vasprun.xml"))
                bands.get_efermi(vasprun=os.path.join(args.directory, "vasprun.xml" if args.efermi.lower() == "nscf" else "vasprun.xml.scf"))
                bands.export(directory=args.directory, engine=args.engine, bandrange=args.bandrange, xrange=args.xrange, yrange=args.yrange)
                if bands.magnetic_status == "soc-ispin-1" or bands.magnetic_status == "soc-ispin-2":
                    # actually soc-ispin-2 never exists
                    print("=====================================================================\n")
                    print("                             postflow\n")
                    print("---------------------------------------------------------------------\n")
                    print("Note:\n")
                    print("even when you set soc and ISPIN=2 at the same time in INCAR\n")
                    print("you will only find band-structure-soc-ispin-1-xxx in post-processing\n")
                    print("because VASP will actually reset ISPIN to 1 when considering soc\n")
                bands.print_gap()
                bands.print_effective_mass()
            elif args.static == "dos":
                from pymatflow.vasp.post.pdos import post_pdos
                pdos = post_pdos()
                pdos.get_vasprun(os.path.join(args.directory, "vasprun.xml"))
                pdos.get_efermi(vasprun=os.path.join(args.directory, "vasprun.xml" if args.efermi.lower() == "nscf" else "vasprun.xml.scf"))            
                pdos.export(directory=args.directory, engine=args.engine, plotrange=args.plotrange)
                
        elif args.runtype == 1:
            # optimization
            from pymatflow.vasp.post.opt import opt_out
            os.chdir(args.directory)
            opt = opt_out()
            opt.get_info(outcar="OUTCAR", poscar="POSCAR")
            os.chdir("../")
            opt.export(directory=args.directory)
        elif args.runtype == 2:
            # cubic cell
            os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
        elif args.runtype == 3:
            # hexagonal cell
            os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
        elif args.runtype == 4:
            # tetragonal cell
            os.system("cd %s; bash get_energy.sh; rm get_energy.sh; cd ../../" % os.path.join(args.directory, "post-processing"))
        elif args.runtype == 5:
            # neb
            pass
        elif args.runtype == 6:
            # vasp phonon
            from pymatflow.vasp.post.phonon import phonon_post
            task = phonon_post()
            task.supercell_n = args.supercell_n
            task.get_kpath(get_kpath(args.kpath_manual, args.kpath_file))
            task.get_xyz(xyzfile)
            task.export(directory=args.directory)
        elif args.runtype == 7:
            # phonopy
            from pymatflow.vasp.post.phonopy import phonopy_post
            task = phonopy_post()
            task.supercell_n = args.supercell_n
            task.get_kpath(get_kpath(args.kpath_manual, args.kpath_file))
            task.get_xyz(xyzfile)
            task.export(directory=args.directory)
    # --------------------------------------------------------------------------



if __name__ == "__main__":
    main()
