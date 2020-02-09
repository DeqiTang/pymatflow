#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from pymatflow.qe.tddfpt import tddfpt_run
from pymatflow.remote.server import server_handle

"""
usage:
"""

lr_input_tl = {}
lr_control_tl = {}
lr_input_ts = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-static",
            help="directory of the calculation")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")


    # ------------------------------------------------------
    #                 lr_input
    # ------------------------------------------------------
    parser.add_argument("--prefix", type=str, default="pwscf",
            help="prefix used in pw.x")

    parser.add_argument("--outdir", type=str, default="./tmp",
            help="outdir used in pw.x")

    # ----------------------------------------------------------------
    #                         lr_control
    # ----------------------------------------------------------------
    parser.add_argument("--itermax", type=int, default=1500,
            help="itermax in lr_control in turbo_lanczos.x calc")

    parser.add_argument("--ipol", type=int, default=1,
            help="ipol in lr_control in turbo_lanczos.x calc")



    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")
    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")
    parser.add_argument("--jobname", type=str, default="turbo_lanczos",
            help="jobname on the pbs server")
    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")
    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

    lr_input_tl["prefix"] = args.prefix
    lr_input_tl["outdir"] = args.outdir
    lr_input_ts["prefix"] = args.prefix
    lr_input_ts["outdir"] = args.outdir

    lr_control_tl["itermax"] =args.itermax
    lr_control_tl["ipol"] = args.ipol


    task = tddfpt_run()
    #task.get_xyz(args.file)
    task.set_turbo_lanczos(lr_input=lr_input_tl, lr_control=lr_control_tl)
    task.set_turbo_spectrum(lr_input=lr_input_ts)
    task.turbo_lanczos(directory=args.directory, runopt=args.runopt, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)

    server_handle(auto=args.auto, directory=args.directory, jobfilebase="turbo-lanczos", server=args.server)
