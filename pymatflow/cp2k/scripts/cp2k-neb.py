#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.cp2k.neb import neb_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync

"""
usage:
    cp2k-scf.py -f xxx.xyz
"""


force_eval = {}
motion = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of the calculation", type=str, default="tmp-cp2k-neb")

    parser.add_argument("--runopt", type=str, default="genrun",
            help="runopt: 'gen', 'run', 'genrun'")

    parser.add_argument("-f", "--file", help="the xyz file name", type=str)

    parser.add_argument("-p", "--printout-option", nargs="+", type=int,
            default=[],
            help="Properties printout option(0, 1, 2 implemented now), you can also activate multiple prinout-option at the same time")

    
    # ------------------------------------------------------------------
    #                    force_eval/dft related parameters
    # ------------------------------------------------------------------
    
    parser.add_argument("--qs-method", type=str, default="gpw",
            choices=["am1", "dftb", "gapw", "gapw_xc", "gpw", "lrigpw", "mndo", "mndod", 
                "ofgpw", "pdg", "pm3", "pm6", "pm6-fm", "pnnl", "rigpw", "rm1"],
            help="dft-qs-method: specify the electronic structure method")

    parser.add_argument("--eps-scf", type=float, default=1.0e-6,
            help="dft-scf-eps_scf")

    parser.add_argument("--xc-functional", type=str, default="pbe",
            help="dft-xc-xc_functional: LYP, PADE, PBE, PW92, TPSS, XGGA, XWPBE, etc.")

    parser.add_argument("--cutoff", type=int, default=100,
            help="CUTOFF, default value: 100 Ry")

    parser.add_argument("--rel-cutoff", type=int, default=60,
            help="REL_CUTOFF, default value: 60 Ry")

    parser.add_argument("-k", "--kpoints-scheme", type=str,
            default="GAMMA",
            help="DFT-KPOINTS-SCHEME(str): can be NONE, GAMMA, MONKHORST-PACK, MACDONALD, GENERAL. when you set MONKHORST-PACK, you should also add the three integers like 'monkhorst-pack 3 3 3'")
    
    parser.add_argument("--diag", type=str, default="TRUE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing tranditional diagonalization for SCF")

    parser.add_argument("--ot", type=str, default="FALSE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing orbital transformation for SCF")

    parser.add_argument("--alpha", type=float, default=0.4,
            help="DFT-SCF-MIXING-ALPHA")

    parser.add_argument("--smear", type=str, default="FALSE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="switch on or off smearing for occupation")
    
    parser.add_argument("--smear-method", type=str, default="FERMI_DIRAC",
            help="smearing type: FERMI_DIRAC, ENERGY_WINDOW")

    parser.add_argument("--added-mos", type=int, default=0,
            help="ADDED_MOS for SCF")


    parser.add_argument("--electronic-temp", type=float, default=300,
            help="ELECTRON_TEMPERATURE for FERMI_DIRAC SMEAR")

    parser.add_argument("--window-size", type=float, default=0,
            help="Size of the energy window centred at the Fermi level for ENERGY_WINDOW type smearing")

    parser.add_argument("--ls-scf", type=str, default="false",
            #choices=["true", "false", "true", "false"],
            help="dft-ls_scf: use linear scaling scf method")

    # vdw correction related
    parser.add_argument("--vdw", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="whether to use VDW correction")

    parser.add_argument("--vdw-potential-type", type=str, default="PAIR_POTENTIAL",
            choices=["PAIR_POTENTIAL", "NON_LOCAL", "NONE"],
            help="DFT-XC-VDW_POTENTIAL-POTENTIAL_TYPE: PAIR_POTENTIAL, NON_LOCAL")

    parser.add_argument("--pair-type", type=str, default="DFTD3",
            choices=["DFTD2", "DFTD3", "DFTD3(BJ)"],
            help="VDW PAIR_POTENTIAL type: DFTD2, DFTD3, DFTD3(BJ)")

    parser.add_argument("--r-cutoff", type=float, default=1.05835442E+001,
            help="DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL: Range of potential. The cutoff will be 2 times this value")


    # --------------------------------
    # MOTION/BAND
    # --------------------------------
    parser.add_argument("--images", nargs="+", type=str,
            help="specify the image xyz file(--images first.xyz imtermediate-1.xyz intermediate-2.xyz ... last.xyz)")

    parser.add_argument("--band-type", type=str, default="CI-NEB",
            help="specify the type of band calculation")

    parser.add_argument("--number-of-replica", type=int, default=5,
            help="number of replicas")

    parser.add_argument("--k-spring", type=float, default=2.0e-2,
            help="value of the spring constant")

    parser.add_argument("--align-frames", type=str, default="TRUE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="Enables the alignment of the frames at the beginning of a BAND calculation. This keyword does not affect the rotation of the replicas during a BAND calculation.")
    
    parser.add_argument("--rotate-frames", type=str, default="TRUE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="Compute at each BAND step the RMSD and rotate the frames in order to minimize it.")

    # for server
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.emuhelper/server.conf")
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()

    force_eval["dft-ls_scf"] = args.ls_scf
    force_eval["dft-qs-method"] = args.qs_method
    force_eval["dft-mgrid-cutoff"] = args.cutoff
    force_eval["dft-mgrid-rel_cutoff"] = args.rel_cutoff
    force_eval["dft-xc-xc_functional"] = args.xc_functional
    force_eval["dft-scf-eps_scf"] = args.eps_scf
    force_eval["dft-scf-added_mos"] = args.added_mos
    force_eval["dft-scf-smear"] = args.smear
    force_eval["dft-scf-smear-method"] = args.smear_method
    force_eval["dft-scf-smear-electronic_temperature"] = args.electronic_temp
    force_eval["dft-scf-smear-window_size"] = args.window_size
    force_eval["dft-scf-diagonalization"] = args.diag
    force_eval["dft-scf-ot"] = args.ot
    force_eval["dft-scf-mixing-alpha"] = args.alpha
    force_eval["dft-kpoints-scheme"] = args.kpoints_scheme

    force_eval["dft-xc-vdw_potential"] = args.vdw
    force_eval["dft-xc-vdw_potential-potential_type"] = args.vdw_potential_type
    force_eval["dft-xc-vdw_potential-pair_potential-type"] = args.pair_type
    force_eval["dft-xc-vdw_potential-pair-potential-r_cutoff"] = args.r_cutoff

    motion["BAND-BAND_TYPE"] = args.band_type
    motion["BAND-NUMBER_OF_REPLICA"] = args.number_of_replica
    motion["BAND-ALIGN_FRAMES"] = args.align_frames
    motion["BAND-ROTATE-FRAMES"] = args.rotate_frames
    motion["BAND-K_SPRING"] = args.k_spring


    task = neb_run(images=args.images)
    task.neb(directory=args.directory, runopt=args.runopt, force_eval=force_eval, motion=motion)

    # server handle
    if args.auto == 0:
        pass
    elif args.auto == 1:
        mover = rsync()
        mover.get_info(os.path.join(os.path.expanduser("~"), ".emuhelper/server.conf"))
        mover.copy_default(source=os.path.abspath(args.directory))
    elif args.auto == 2:
        mover = rsync()
        mover.get_info(os.path.join(os.path.expanduser("~"), ".emuhelper/server.conf"))
        mover.copy_default(source=os.path.abspath(args.directory))
        ctl = ssh()
        ctl.get_info(os.path.join(os.path.expanduser('~'), ".emuhelper/server.conf"))
        ctl.login()
        ctl.submit(workdir=args.directory, jobfile="neb.inp.sub")
