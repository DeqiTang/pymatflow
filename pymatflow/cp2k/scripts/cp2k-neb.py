#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.cp2k.neb import neb_run
from pymatflow.remote.server import server_handle

"""
usage:
    cp2k-scf.py -f xxx.xyz
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-cp2k-neb",
            help="directory to do the Nudged Elastic Band calculation")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("-f", "--file",type=str,
            help="the xyz strcture file with second line specifying cell parameters")

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

    parser.add_argument("--smear-method", type=str, default="FERMI_DIRAC",
            help="smearing type: FERMI_DIRAC, ENERGY_WINDOW")

    parser.add_argument("--added-mos", type=int, default=0,
            help="ADDED_MOS for SCF")

    parser.add_argument("--electronic-temp", type=float, default=300,
            help="ELECTRON_TEMPERATURE for FERMI_DIRAC SMEAR")

    parser.add_argument("--window-size", type=float, default=0,
            help="Size of the energy window centred at the Fermi level for ENERGY_WINDOW type smearing")

    parser.add_argument("--ls-scf", type=str, default="false",
            choices=["true", "false", "true", "false"],
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


    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="cp2-neb",
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

    params["FORCE_EVAL-dft-xc-vdw_potential"] = args.vdw
    params["FORCE_EVAL-dft-xc-vdw_potential-potential_type"] = args.vdw_potential_type
    params["FORCE_EVAL-dft-xc-vdw_potential-pair_potential-type"] = args.pair_type
    params["FORCE_EVAL-dft-xc-vdw_potential-pair-potential-r_cutoff"] = args.r_cutoff

    params["MOTION-BAND-BAND_TYPE"] = args.band_type
    params["MOTION-BAND-NUMBER_OF_REPLICA"] = args.number_of_replica
    params["MOTION-BAND-ALIGN_FRAMES"] = args.align_frames
    params["MOTION-BAND-ROTATE-FRAMES"] = args.rotate_frames
    params["MOTION-BAND-K_SPRING"] = args.k_spring


    task = neb_run()
    task.get_images(images=args.images)
    task.set_params(params=params)
    task.neb(directory=args.directory, runopt=args.runopt, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)

    server_handle(auto=args.auto, directory=args.directory, jobfilebase="neb", server=args.server)
