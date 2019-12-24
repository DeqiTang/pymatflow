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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()    
    parser.add_argument("-d", "--directory", help="directory for the static running", type=str, default="tmp-qe-static")
    parser.add_argument("--mpi", help="MPI commadn", type=str, default="")
    parser.add_argument("-f", "--file", help="the xyz file", type=str)
    parser.add_argument("--runopt", help="gen, run, or genrun", type=str, default="genrun")
   
    # --------------------------------------------------------------
    # for dynmat
    # --------------------------------------------------------------
    parser.add_argument("--q", type=float, nargs="+",
            default=[0, 0, 0],
            help="calculate LO modes along the direction q.")

    parser.add_argument("--asr", type=str, default="simple",
            choices=["no", "simple", "crystal", "one-dim", "zero-dim"],
            help="the type of Acoustic Sum Rule imposed")

    # --------------------------------------------------------------
    # for server
    # --------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.emuhelper/server.conf")
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file


    #task = static_run(xyzfile)
    task = dfpt_run(xyzfile)
    task.dynmat(directory=args.directory, mpi=args.mpi, runopt=args.runopt, asr=args.asr, qi=args.q)

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
        ctl.submit(workdir=args.directory, jobfile="dynmat.in.sub")
        # cannot submit the following job before finishing the previous one
        #ctl.submit(wordir=args.directory, jobfile="q2r.in.sub")
        #ctl.submit(wordir=args.directory, jobfile="matdyn.in.sub")
        #ctl.submit(wordir=args.directory, jobfile="plotband.in.sub")

