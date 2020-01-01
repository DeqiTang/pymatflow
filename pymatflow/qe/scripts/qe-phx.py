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
    parser.add_argument("--tr2-ph", type=float, default=1.0e-14,
            help="threshold for self-consistency.")

    parser.add_argument("--nq", type=int, nargs="+",
            default=[0, 0, 0],
            help="set value of nq1 nq2 nq3.")

    parser.add_argument("--epsil", type=str, default=None,
            choices=[".true.", ".false."],
            help="set epsil in inputph")
    
    parser.add_argument("--lraman", type=str, default=None,
            choices=["true", "false"],
            help="set lraman, can be 'true' or 'false' only. default is None which means 'false' in real world.")

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
    inputph["lraman"] = args.lraman
    inputph["epsil"] = args.epsil
    inputph["nq1"] = args.nq[0]
    inputph["nq2"] = args.nq[1]
    inputph["nq3"] = args.nq[2]


    task = dfpt_run()
    task.get_xyz(xyzfile)
    task.set_inputph(inputph=inputph)
    task.phx(directory=args.directory, mpi=args.mpi, runopt=args.runopt)

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
        ctl.submit(workdir=args.directory, jobfile="phx.sub")
        # cannot submit the following job before finishing the previous one
        #ctl.submit(wordir=args.directory, jobfile="q2r.in.sub")
        #ctl.submit(wordir=args.directory, jobfile="matdyn.in.sub")
        #ctl.submit(wordir=args.directory, jobfile="plotband.in.sub")

