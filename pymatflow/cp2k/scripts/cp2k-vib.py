#!/usr/bin/env python
# _*_ coding: utf-8 _*_


import os
import argparse

from pymatflow.cp2k.vib import vib_run
from pymatflow.remote.server import server_handle
"""
usage:
    cp2k-vib.py -f xxx.xyz
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-cp2k-vib",
            help="directory to do the vib calculation")

    parser.add_argument("-f", "--file", type=str,
            help="the xyz structure file with second line specifying cell parameter")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")


    # ---------------------------------------------------------------------------
    #                       FORCE_EVAL realated parameters
    # ---------------------------------------------------------------------------

    parser.add_argument("--ls-scf", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="use linear scaling scf method")

    parser.add_argument("--qs-method", type=str, default="GPW",
            choices=["AM1", "DFTB", "GAPW", "GAPW_XC", "GPW", "LRIGPW", "MNDO", "MNDOD",
                "OFGPW", "PDG", "PM3", "PM6", "PM6-FM", "PNNL", "RIGPW", "RM1"],
            help="specify the electronic structure method")

    parser.add_argument("--eps-scf", type=float, default=1.0e-6,
            help="DFT-SCF-EPS_SCF")

    parser.add_argument("--xc-functional", type=str, default="PBE",
            help="XC_FUNCTIONAL type")

    parser.add_argument("--cutoff", type=int, default=100,
            help="CUTOFF, default value: 100 Ry")

    parser.add_argument("--rel-cutoff", type=int, default=60,
            help="REL_CUTOFF, default value: 60 Ry")

    parser.add_argument("-k", "--kpoints-scheme", type=str,
            default="GAMMA",
            help="DFT-KPOINTS-SCHEME(str): can be NONE, GAMMA, MONKHORST-PACK, MACDONALD, GENERAL. when you set MONKHORST-PACK, you should also add the three integers like 'monkhorst-pack 3 3 3'")

    parser.add_argument("--diag", type=str, default="TRUE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing tranditional diagonalization for SCF")

    parser.add_argument("--ot", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing orbital transformation for SCF")

    parser.add_argument("--alpha", type=float, default=0.4,
            help="DFT-SCF-MIXING-ALPHA")

    parser.add_argument("--smear", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="switch on or off smearing for occupation")

    parser.add_argument("--added-mos", type=int, default=0,
            help="ADDED_MOS for SCF")

    parser.add_argument("--smear-method", type=str, default="FERMI_DIRAC",
            help="smearing type: FERMI_DIRAC, ENERGY_WINDOW")

    parser.add_argument("--electronic-temp", type=float, default=300,
            help="ELECTRON_TEMPERATURE for FERMI_DIRAC SMEAR")

    parser.add_argument("--window-size", type=float, default=0,
            help="Size of the energy window centred at the Fermi level for ENERGY_WINDOW type smearing")



    # ---------------------------------------------------------------
    #             vibrational_analysis related parameters
    # ---------------------------------------------------------------

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

    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="geo-opt",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

    params = {}

    params["FORCE_EVAL-DFT-LS_SCF"] = args.ls_scf
    params["FORCE_EVAL-DFT-QS-METHOD"] = args.qs_method
    params["FORCE_EVAL-DFT-MGRID-CUTOFF"] = args.cutoff
    params["FORCE_EVAL-DFT-MGRID-REL_CUTOFF"] = args.rel_cutoff
    params["FORCE_EVAL-DFT-XC-XC_FUNCTIONAL"] = args.xc_functional
    params["FORCE_EVAL-DFT-SCF-EPS_SCF"] = args.eps_scf
    params["FORCE_EVAL-DFT-SCF-ADDED_MOS"] = args.added_mos
    params["FORCE_EVAL-DFT-SCF-SMEAR"] = args.smear
    params["FORCE_EVAL-DFT-SCF-SMEAR-METHOD"] = args.smear_method
    params["FORCE_EVAL-DFT-SCF-SMEAR-ELECTRONIC_TEMPERATURE"] = args.electronic_temp
    params["FORCE_EVAL-DFT-SCF-SMEAR-WINDOW_SIZE"] = args.window_size
    params["FORCE_EVAL-DFT-SCF-DIAGONALIZATION"] = args.diag
    params["FORCE_EVAL-DFT-SCF-OT"] = args.ot
    params["FORCE_EVAL-DFT-SCF-MIXING-ALPHA"] = args.alpha
    params["FORCE_EVAL-DFT-KPOINTS-SCHEME"] = args.kpoints_scheme

    params["VIBRATIONAL_ANALYSIS-DX"] = args.dx
    params["VIBRATIONAL_ANALYSIS-FULLY_PERIODIC"] = args.fully_periodic
    params["VIBRATIONAL_ANALYSIS-INTENSITIES"] = args.intensities
    params["VIBRATIONAL_ANALYSIS-TC_PRESSURE"] = args.tc_pressure
    params["VIBRATIONAL_ANALYSIS-TC_TEMPERATURE"] = args.tc_temperature
    params["VIBRATIONAL_ANALYSIS-THERMOCHEMISTRY"] = args.thermochemistry

    task = vib_run()
    task.get_xyz(args.file)
    task.set_params(params=params)
    task.vib(directory=args.directory, runopt=args.runopt, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)

    server_handle(auto=args.auto, directory=args.directory, jobfilebase="vib", server=args.server)
