#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from pymatflow.qe.tddfpt import tddfpt_run


"""
usage:
"""

lr_input_td = {}
lr_dav_td = {}
lr_input_ts = {}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-static",
            help="directory of the calculation")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    # ---------------------------------------------------------------
    #                       lr_input
    # ---------------------------------------------------------------
    parser.add_argument("--prefix", type=str, default="pwscf",
            help="prefix used in pw.x")

    parser.add_argument("--outdir", type=str, default="./tmp",
            help="outdir used in pw.x")

    # -------------------------------------------------------------
    #                         lr_dav
    # -------------------------------------------------------------
    parser.add_argument("--if-dft-spectrum", type=str, default="False",
            choices=["True", "False"],
            help="if_dft_spectrum in lr_dav in turbo_davidson.x calc")

    parser.add_argument("--num-eign", type=int, default=15,
            help="num_eign in lr_dav in turbo_davidson.x calc")

    parser.add_argument("--num-init", type=int, default=30,
            help="num_init in lr_dav in turbo_davidson.x calc")

    parser.add_argument("--num-basis-max", type=int, default=90,
            help="num_basis_max in lr_dav in turbo_davidson.x calc")

    parser.add_argument("--residue-conv-thr", type=float, default=1.0e-6,
            help="residue_conv_thr in lr_dav in turbo_davidson.x calc")

    parser.add_argument("--start", type=float, default=0.0,
            help="start in residue_conv_thr in turbo_davidson.x calc")

    parser.add_argument("--finish", type=float, default=1.0,
            help="finish in lr_dav in turbo_davidson.x calc")

    parser.add_argument("--step", type=float, default=0.01,
            help="step in lr_dav in turbo_davidson.x calc")

    parser.add_argument("--broadening", type=float, default=0.04,
            help="broadenging in turbo_davidson.x calc")

    parser.add_argument("--reference", type=float, default=0.04,
            help="reference in turbo_davidson.x calc")

    # -----------------------------------------------------------------
    #                       run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="turbo_davidson",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

    lr_input_td["prefix"] = args.prefix
    lr_input_td["outdir"] = args.outdir
    lr_input_ts["prefix"] = args.prefix
    lr_input_ts["outdir"] = args.outdir

    lr_dav_td["if_dft_spectrum"] = args.if_dft_spectrum
    lr_dav_td["num_eign"] = args.num_eign
    lr_dav_td["num_init"] = args.num_init
    lr_dav_td["num_basis_max"] = args.num_basis_max
    lr_dav_td["residue_conv_thr"] = args.residue_conv_thr
    lr_dav_td["start"] = args.start
    lr_dav_td["finish"] = args.finish
    lr_dav_td["step"] = args.step
    lr_dav_td["broadening"] = args.broadening
    lr_dav_td["reference"] = args.reference

    task = tddfpt_run()
    #task.get_xyz(args.file)
    task.set_turbo_davidson(lr_input=lr_input_td, lr_dav=lr_dav_td)
    task.set_turbo_spectrum(lr_input=lr_input_ts)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.turbo_davidson(directory=args.directory, runopt=args.runopt, auto=args.auto)
