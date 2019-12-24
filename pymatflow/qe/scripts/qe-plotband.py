#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

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

    task = dfpt_run(xyzfile)
    task.plotband(directory=args.directory, mpi=args.mpi, runopt=args.runopt, freq_min=args.freq_min, freq_max=args.freq_max, efermi=args.efermi, freq_step=args.freq_step, freq_reference=args.freq_reference)

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
        ctl.submit(workdir=args.directory, jobfile="plotband.in.sub")
        # cannot submit the following job before finishing the previous one
        #ctl.submit(wordir=args.directory, jobfile="q2r.in.sub")
        #ctl.submit(wordir=args.directory, jobfile="matdyn.in.sub")
        #ctl.submit(wordir=args.directory, jobfile="plotband.in.sub")

