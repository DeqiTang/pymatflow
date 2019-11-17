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
            choicies = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
            help="Properties printout option, you can also activate multiple prinout-option at the same time. 
            1: printout pdos
            2: printout band
            3: printout electron densities
            4: printout electron local function(ELF)
            5: printout molecular orbitals
            6: printout molecular orbital cube files
            7: printout mulliken populaltion analysis
            8: printout cubes for generation of STM images
            9: printout cube file with total density(electrons+atomic core)
           10: printout v_hartree_cube
           11: printout v_xc_cube
           12: printout xray_diffraction_spectrum
           13: request a RESP fit of charges.")

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

    parser.add_argument("--smear", type=str, default="FALSE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="switch on or off smearing for occupation")

    parser.add_argument("--added-mos", help="ADDED_MOS for SCF", type=int, default=0)

    parser.add_argument("--smear-method", help="smearing type: FERMI_DIRAC, ENERGY_WINDOW", type=str, default="FERMI_DIRAC")

    parser.add_argument("--electronic-temp", help="ELECTRON_TEMPERATURE for FERMI_DIRAC SMEAR", type=float, default=300)

    parser.add_argument("--window-size", help="Size of the energy window centred at the Fermi level for ENERGY_WINDOW type smearing", type=float, default=0)

    # vdw correction related
    parser.add_argument("--vdw", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="whether to use VDW correction")

    parser.add_argument("--vdw-potential-type", type=str, default="PAIR_POTENTIAL",
            choices=["PAIR_POTENTIAL", "NON_LOCAL", "NONE"],
            help="type of VDW POTENTIAL: PAIR_POTENTIAL, NON_LOCAL")

    parser.add_argument("--pair-type", type=str, default="DFTD3",
            choices=["DFTD2", "DFTD3", "DFTD3(BJ)"],
            help="VDW PAIR_POTENTIAL type: DFTD2, DFTD3, DFTD3(BJ)")

    parser.add_argument("--r-cutoff", type=float, default=1.05835442E+001,
            help="Range of potential. The cutoff will be 2 times this value")

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

    force_eval["DFT-XC-VDW_POTENTIAL"] = args.vdw
    force_eval["DFT-XC-VDW_POTENTIAL-POTENTIAL_TYPE"] = args.vdw_potential_type
    force_eval["DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-TYPE"] = args.pair_type
    force_eval["DFT-XC-VDW_POTENTIAL-PAIR-POTENTIAL-R_CUTOFF"] = args.r_cutoff

    task = static_run(xyzfile)
    task.scf(directory=directory, runopt=args.runopt, force_eval=force_eval, printout_option=args.printout_option)

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
