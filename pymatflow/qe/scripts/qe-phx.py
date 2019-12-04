#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

#from pymatflow.qe.static import static_run
from pymatflow.qe.dfpt import dfpt_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync

"""
usage:
"""

inputph = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-static",
            help="directory for the static running(also where dpft will be going)")

    parser.add_argument("--mpi", type=str, default="",
            help="MPI commadn")

    parser.add_argument("-f", "--file", type=str,
            help="the xyz file containing the structure")

    parser.add_argument("--runopt", type=str, default="genrun",
            help="gen, run, or genrun")
   
    # --------------------------------------------------------------
    # for phx
    # --------------------------------------------------------------
    parser.add_argument("--qpoints-option", type=str, default="gamma",
            choices=["gamma", "qmesh"],
            help="qpoints option for ph.x")

    parser.add_argument("--qpoints", type=int, nargs="+",
            default=[2, 2, 2],
            help="qpoints mesh for ph.x")

    parser.add_argument("--polar", type=str, default="no",
            choices=["yes", "no"],
            help="whether it is polar materials, if yes, will set epsil = .true. and calculate and store the dielectric tensor and effective charges.")

    parser.add_argument("--tr2-ph", type=float, default=1.0e-14,
            help="threshold for self-consistency.")

    # ------------------------------
    # for server
    # ------------------------------
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.emuhelper/server.conf")
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file

    inputph["tr2_ph"] = args.tr2_ph

    task = dfpt_run(xyzfile)
    task.phx(directory=args.directory, mpi=args.mpi, runopt=args.runopt, qpoints_option=args.qpoints_option, qpoints=args.qpoints, polar=args.polar, inputph=inputph)

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
        ctl.submit(workdir=args.directory, jobfile="phx.in.sub")
        # cannot submit the following job before finishing the previous one
        #ctl.submit(wordir=args.directory, jobfile="q2r.in.sub")
        #ctl.submit(wordir=args.directory, jobfile="matdyn.in.sub")
        #ctl.submit(wordir=args.directory, jobfile="plotband.in.sub")

