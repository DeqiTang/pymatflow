#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from pymatflow.qe.static import StaticRun


"""
usage qe-converge-degauss.py -f xxx.xyz --range degauss_min degauss_max step
"""

control_params = {}
# do not set occupation, smearing and degauss through system_params
system_params = {}
electrons_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--xyz", help="the xyz file name", type=str)

    parser.add_argument("--range", help="degauss_min degauss_max step", nargs='+', type=float)

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    parser.add_argument("--kpoints-option", type=str, default="automatic",
            choices=["automatic", "gamma", "crystal_b"],
            help="Kpoints generation scheme option for the SCF or non-SCF calculation")

    parser.add_argument("--kpoints-mp", type=int, nargs="+",
            default=[1, 1, 1, 0, 0, 0],
            help="Monkhorst-Pack kpoint grid, in format like --kpoints-mp 1 1 1 0 0 0")

    parser.add_argument("--ecutwfc", help="better a previously converged ecutwfc", type=int, default=100)

    parser.add_argument("--ecutrho", type=int, default=None,
            help="Kinetic energy cutoff for charge density and potential in unit of Rydberg, default value: None")

    # -----------------------------------------------------------------
    #                       run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="converge-degauss",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
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
    xyzfile = args.xyz
    system_params["ecutwfc"] = args.ecutwfc
    system_params["ecutrho"] = args.ecutrho

    task = StaticRun()
    task.get_xyz(xyzfile)
    task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=kpoints_mp)
    task.set_params(control=control_params, system=system_params, electrons=electrons_params)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
    task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
    task.converge_degauss(round(args.range[0], 6), round(args.range[1], 6), round(args.range[2], 6), runopt=args.runopt, auto=args.auto)
