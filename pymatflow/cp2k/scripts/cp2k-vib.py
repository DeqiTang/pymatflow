#!/usr/bin/env python
# _*_ coding: utf-8 _*_


import os
import argparse

from pymatflow.cp2k.vib import vib_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync

"""
usage:
    cp2k-vib.py -f xxx.xyz
"""


force_eval = {}
vib = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of the vib calculation", type=str, default="tmp-cp2k-vib")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)

    parser.add_argument("--runopt", type=str, default="genrun", 
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--ls-scf", type=str, default="FALSE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="use linear scaling scf method")
    
    parser.add_argument("--qs-method", type=str, default="GPW",
            choices=["AM1", "DFTB", "GAPW", "GAPW_XC", "GPW", "LRIGPW", "MNDO", "MNDOD", 
                "OFGPW", "PDG", "PM3", "PM6", "PM6-FM", "PNNL", "RIGPW", "RM1"],
            help="specify the electronic structure method")

    parser.add_argument("--eps-scf", type=float, default=1.0e-5,
            help="DFT-SCF-EPS_SCF")

    parser.add_argument("--xc-functional", type=str, default="PBE",
            help="XC_FUNCTIONAL type")

    parser.add_argument("--cutoff", type=int, default=100,
            help="CUTOFF, default value: 100 Ry")

    parser.add_argument("--rel-cutoff", type=int, default=60,
            help="REL_CUTOFF, default value: 60 Ry")

    parser.add_argument("-k", "--kpoints", help="set kpoints like '3 3 3 0 0 0'", type=str, default="3 3 3 0 0 0")
    
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

    parser.add_argument("--added-mos", help="ADDED_MOS for SCF", type=int, default=0)

    parser.add_argument("--smear-method", help="smearing type: FERMI_DIRAC, ENERGY_WINDOW", type=str, default="FERMI_DIRAC")

    parser.add_argument("--electronic-temp", help="ELECTRON_TEMPERATURE for FERMI_DIRAC SMEAR", type=float, default=300)

    parser.add_argument("--window-size", help="Size of the energy window centred at the Fermi level for ENERGY_WINDOW type smearing", type=float, default=0)

    # vibrational_analysis
    parser.add_argument("--dx", type=float, default=1.0e-2,
            help="specify the increment to be used to construct the HESSIAN with finite difference method")

    parser.add_argument("--fully-periodic", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="avoids to clean rotations from the Hessian matrix")

    parser.add_argument("--intensities", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="Calculation of the IR-Intensities. Calculation of dipoles has to be specified explicitly"
            )

    parser.add_argument("--tc-pressure", type=float, default=1.01325000E+005,
            help="Pressure for the calculation of the thermochemical data in unit of [Pa]")

    parser.add_argument("--tc-temperature", type=float, default=2.73150000E+002,
            help="Temperature for the calculation of the thermochemical data in unit of [K]")

    parser.add_argument("--thermochemistry", type=str, default="FALSE",
            help="Calculation of the thermochemical data. Valid for molecules in the gas phase.")

    # for server
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.emuhelper/server.conf")
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()
    directory = args.directory
    xyzfile = args.file
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]

    force_eval["DFT-LS_SCF"] = args.ls_scf
    force_eval["DFT-QS-METHOD"] = args.qs_method
    force_eval["DFT-MGRID-CUTOFF"] = args.cutoff
    force_eval["DFT-MGRID-REL_CUTOFF"] = args.rel_cutoff
    force_eval["DFT-XC-XC_FUNCTIONAL"] = args.xc_functional
    force_eval["DFT-SCF-EPS_SCF"] = args.eps_scf
    force_eval["DFT-SCF-ADDED_MOS"] = args.added_mos
    force_eval["DFT-SCF-SMEAR"] = args.smear
    force_eval["DFT-SCF-SMEAR-METHOD"] = args.smear_method
    force_eval["DFT-SCF-SMEAR-ELECTRONIC_TEMPERATURE"] = args.electronic_temp
    force_eval["DFT-SCF-SMEAR-WINDOW_SIZE"] = args.window_size
    force_eval["DFT-SCF-DIAGONALIZATION"] = args.diag
    force_eval["DFT-SCF-OT"] = args.ot
    force_eval["DFT-SCF-MIXING-ALPHA"] = args.alpha

    vib["DX"] = args.dx
    vib["FULLY_PERIODIC"] = args.fully_periodic
    vib["INTENSITIES"] = args.intensities
    vib["TC_PRESSURE"] = args.tc_pressure
    vib["TC_TEMPERATURE"] = args.tc_temperature
    vib["THERMOCHEMISTRY"] = args.thermochemistry

    task = vib_run(xyzfile)
    task.vib(directory=directory, runopt=args.runopt, force_eval=force_eval, vibrational=vib)

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
        ctl.submit(workdir=args.directory, jobfile="vib.inp.sub")
