#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from pymatflow.qe.tddfpt import tddfpt_run


"""
usage:
"""

inputpp = {}
energy_grid = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-static",
            help="directory of the calculation")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    # ----------------------------------------------------------
    #                        INPUTPP
    # ----------------------------------------------------------
    parser.add_argument("--calculation", type=str, default="eps",
            help="calculation in inputpp for epsilon.x")

    parser.add_argument("--prefix", type=str, default="pwscf",
            help="prefix used in pw.x")

    parser.add_argument("--outdir", type=str, default="./tmp",
            help="outdir used in pw.x")

    # ----------------------------------------------------------
    #                            ENERGY_GRID
    # ----------------------------------------------------------
    parser.add_argument("--smeartype", type=str, default="gaussian",
            help="smeartype in ENERGY_GRID in epsilon.x calc")

    parser.add_argument("--intersmear", type=float, default=0.1,
            help="intersmear in ENERGY_GRID in epsilon.x calc")

    parser.add_argument("--wmin", type=float, default=0.0,
            help="wmin in ENERGY_GRID in epsilon.x calc")

    parser.add_argument("--wmax", type=float, default=15.0,
            help="wmax in ENERGY_GRID in epsilon.x calc")

    parser.add_argument("--nw", type=int, default=1000,)

    # -----------------------------------------------------------------
    #                       run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="qe-epsilon",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

    inputpp["calculation"] = args.calculation
    inputpp["prefix"] = args.prefix
    inputpp["outdir"] = args.outdir

    energy_grid["smeartype"] = args.smeartype
    energy_grid["intersmear"] = args.intersmear
    energy_grid["wmin"] = args.wmin
    energy_grid["wmax"] = args.wmax
    energy_grid["nw"] = args.nw

    task = tddfpt_run()
    #task.get_xyz(args.file)
    task.set_epsilon(inputpp=inputpp, energy_grid=energy_grid)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.epsilon(directory=args.directory, runopt=args.runopt, auto=args.auto)
