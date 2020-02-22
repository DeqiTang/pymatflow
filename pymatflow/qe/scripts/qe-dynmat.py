#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

#from pymatflow.qe.static import static_run
from pymatflow.qe.dfpt import dfpt_run


"""
usage:
"""

dynmat_input = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-static",
            help="directory for the static running")

    parser.add_argument("--mpi", type=str, default="",
            help="MPI commadn")

    parser.add_argument("-f", "--file", type=str,
            help="the xyz file")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    # --------------------------------------------------------------
    # for dynmat
    # --------------------------------------------------------------
    parser.add_argument("--fildyn", type=str, default="phx.dyn",
            help="specify fildyn which contains frequency info.")

    parser.add_argument("--qi", type=float, nargs="+",
            default=[0, 0, 0],
            help="calculate LO modes along the direction q.")

    parser.add_argument("--asr", type=str, default="simple",
            choices=["no", "simple", "crystal", "one-dim", "zero-dim"],
            help="the type of Acoustic Sum Rule imposed")

    # -----------------------------------------------------------------
    #                       run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="qe-dynmat",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

    dynmat_input["fildyn"] = args.fildyn
    dynmat_input["asr"] = args.asr
    dynmat_input["q(1)"] = args.qi[0]
    dynmat_input["q(2)"] = args.qi[1]
    dynmat_input["q(3)"] = args.qi[2]

    task = dfpt_run()
    task.get_xyz(args.file)
    task.set_dynmat(dynmat_input=dynmat_input)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.dynmat(directory=args.directory, mpi=args.mpi, runopt=args.runopt, auto=args.auto)
