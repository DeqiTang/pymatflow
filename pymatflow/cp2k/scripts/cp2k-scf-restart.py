#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.cp2k.static import static_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync

"""
usage:
    cp2k-scf.py -f xxx.xyz
"""


force_eval = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of the calculation", type=str, default="tmp-cp2k-static")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)

    parser.add_argument("--runopt", type=str, default="genrun", 
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("-p", "--printout-option", nargs="+", type=int,
            default=[],
            help="Properties printout option(0, 1, 2 implemented now), you can also activate multiple prinout-option at the same time")
    parser.add_argument("--xc-functional", type=str, default="PBE",
            help="XC_FUNCTIONAL type")
    parser.add_argument("--cutoff", type=int, default=100,
            help="CUTOFF, default value: 100 Ry")
    parser.add_argument("--rel-cutoff", type=int, default=60,
            help="REL_CUTOFF, default value: 60 Ry")

    parser.add_argument("-k", "--kpoints-scheme", type=str,
            default="GAMMA",
            help="DFT-KPOINTS-SCHEME(str): can be NONE, GAMMA, MONKHORST-PACK, MACDONALD, GENERAL. when you set MONKHORST-PACK, you should also add the three integers like 'monkhorst-pack 3 3 3'")

    parser.add_argument("--added-mos", help="ADDED_MOS for SCF", type=int, default=0)
    parser.add_argument("--smear", type=bool, default=False,
            choices=[True, False],
            help="switch on or off smearing for occupation")
    parser.add_argument("--smear-method", help="smearing type: FERMI_DIRAC, ENERGY_WINDOW", type=str, default="FERMI_DIRAC")
    parser.add_argument("--electronic-temp", help="ELECTRON_TEMPERATURE for FERMI_DIRAC SMEAR", type=float, default=300)
    parser.add_argument("--window-size", help="Size of the energy window centred at the Fermi level for ENERGY_WINDOW type smearing", type=float, default=0)

    # for server
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.emuhelper/server.conf")
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()

    force_eval["DFT-MGRID-CUTOFF"] = args.cutoff
    force_eval["DFT-MGRID-REL_CUTOFF"] = args.rel_cutoff
    force_eval["DFT-XC-XC_FUNCTIONAL"] = args.xc_functional
    force_eval["DFT-SCF-ADDED_MOS"] = args.added_mos
    force_eval["DFT-SCF-SMEAR"] = args.smear
    force_eval["DFT-SCF-SMEAR-METHOD"] = args.smear_method
    force_eval["DFT-SCF-SMEAR-ELECTRONIC_TEMPERATURE"] = args.electronic_temp
    force_eval["DFT-SCF-SMEAR-WINDOW_SIZE"] = args.window_size
    force_eval["DFT-KPOINTS-SCHEME"] = args.kpoints_scheme

    task = static_run(args.file)
    task.scf_restart(directory=args.directory, runopt=args.runopt, force_eval=force_eval)

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
        ctl.submit(workdir=args.directory, jobfile="static-scf.inp.sub")
