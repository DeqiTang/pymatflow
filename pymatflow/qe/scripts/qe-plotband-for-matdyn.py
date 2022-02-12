#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.dfpt import dfpt_run

"""
usage:
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-static",
            help="directory for the static running")

    parser.add_argument("--xyz", type=str,
            help="the xyz file")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    # --------------------------------------------------------------
    # for plotband
    # --------------------------------------------------------------
    parser.add_argument("--freq-min", type=float, default=0,
            help="range of frequencies for visualization")

    parser.add_argument("--freq-max", type=float, default=600,
            help="range of frequencies for visualization")

    parser.add_argument("--efermi", type=float, default=0,
            help="fermi energy level(only needed for band structure plot)")

    parser.add_argument("--freq-step", type=float, default=100.0,
            help="freq step")

    parser.add_argument("--freq-reference", type=float, default=0.0,
            help="freq reference")

    # -----------------------------------------------------------------
    #                       run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="plotband-for-matdyn",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            choices=[1], # plotband.x can only run on one node
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    parser.add_argument("--queue", type=str, default=None,
            help="the queue to submit to job, default is not set")

    # llhpc
    parser.add_argument("--partition", type=str, default="free",
            help="choose partition to submit job")

    parser.add_argument("--ntask", type=int, default=24,
            help="choose task number")

    parser.add_argument("--stdout", type=str, default="slurm.out",
            help="set standard out")

    parser.add_argument("--stderr", type=str, default="slurm.err",
            help="set standard err")



    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

    task = dfpt_run()
    task.get_xyz(args.xyz)
    task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)       
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=qrgs.queue)
    task.plotband_for_matdyn(directory=args.directory, runopt=args.runopt, auto=args.auto, freq_min=args.freq_min, freq_max=args.freq_max, efermi=args.efermi, freq_step=args.freq_step, freq_reference=args.freq_reference)
