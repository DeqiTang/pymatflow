#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import argparse

from pymatflow.flow.surface_pes import cp2k_run


"""
usage:
"""

params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-cp2k-pes-opt",
            help="directory of the calculation")

    parser.add_argument("-f", "--file", type=str,
            help="the xyz file name")

    parser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("-p", "--printout-option", nargs="+", type=int,
            default=[],
            help="Properties printout option(0, 1, 2 implemented now), you can also activate multiple prinout-option at the same time")

    parser.add_argument("--ls-scf", type=str, default="FALSE",
            #choices=["TRUE", "FALSE", "true", "false"],
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
            help="CUTOFF, default value: 100 Ry, if you find your SCF hard to converge, you can try increasing CUTOFF")

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

    # motion/geo_opt related
    parser.add_argument("--optimizer", type=str, default="BFGS",
            help="optimization algorithm for geometry optimization: BFGS, CG, LBFGS")

    parser.add_argument("--max-iter", type=int, default=200,
            help="maximum number of geometry optimization steps.")

    parser.add_argument("--type", type=str, default="MINIMIZATION",
            help="specify which kind of geometry optimization to perform: MINIMIZATION(default), TRANSITION_STATE")

    parser.add_argument("--max-dr", type=float, default=3e-3,
            help="Convergence criterion for the maximum geometry change between the current and the last optimizer iteration.")

    parser.add_argument("--max-force", type=float, default=4.50000000E-004,
            help="Convergence criterion for the maximum force component of the current configuration.")

    parser.add_argument("--rms-dr", type=float, default=1.50000000E-003,
            help="Convergence criterion for the root mean square (RMS) geometry change between the current and the last optimizer iteration.")

    parser.add_argument("--rms-force", type=float, default=3.00000000E-004,
            help="Convergence criterion for the root mean square (RMS) force of the current configuration.")


    # ---------------------------------------------------------------
    # for PES
    # ---------------------------------------------------------------
    parser.add_argument("--move-atom", type=int, nargs="+",
            default=[-1],
            help="specify the atoms to move.")

    parser.add_argument("--xrange", type=float, nargs="+",
            default=[1, 3, 0.5],
            help="x range for moving the specified moving atoms.")

    parser.add_argument("--yrange", type=float, nargs="+",
            default=[3, 5, 0.5],
            help="y range for moving the specified moving atoms.")

    parser.add_argument("--zshift", type=float,
            default=0.0,
            help="z shift for the moving atoms, will shift the z of specified moving atoms by value of zshift")

    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")
    parser.add_argument("--jobname", type=str, default="pes-opt",
            help="jobname on the pbs server")
    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")
    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")




    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

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

    params["MOTION-GEO_OPT-MAX_ITER"] = args.max_iter
    params["MOTION-GEO_OPT-OPTIMIZER"] = args.optimizer
    params["MOTION-GEO_OPT-TYPE"] = args.type
    params["MOTION-GEO_OPT-MAX_DR"] = args.max_dr
    params["MOTION-GEO_OPT-MAX_FORCE"] = args.max_force
    params["MOTION-GEO_OPT-RMS_DR"] = args.rms_dr
    params["MOTION-GEO_OPT-RMS_FORCE"] = args.rms_force

    task = cp2k_run()
    task.get_xyz(args.file)
    task.set_params(params=params)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.set_pes( move_atom=args.move_atom, xrange=args.xrange, yrange=args.yrange, zshift=args.zshift)
    task.run(directory=args.directory, runopt=args.runopt, auto=args.auto)
