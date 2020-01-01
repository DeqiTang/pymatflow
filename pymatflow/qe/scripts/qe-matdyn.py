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

matdyn_input = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()    
    parser.add_argument("-d", "--directory", help="directory for the static running", type=str, default="tmp-qe-static")
    parser.add_argument("--mpi", help="MPI commadn", type=str, default="")
    parser.add_argument("-f", "--file", help="the xyz file", type=str)
    parser.add_argument("--runopt", help="gen, run, or genrun", type=str, default="genrun")
   
    # --------------------------------------------------------------
    # for matdyn
    # --------------------------------------------------------------
    parser.add_argument("--asr", type=str, default='simple',
            help="type of sum rule")
    #parser.add_argument("--nqpoints", type=int, default=2,
    #        help="number of qpoints")
    #parser.add_argument("--matdyn-qpoints", type=float, nargs="+",
    #        default= [0.0, 0.0, 0.0, 0.0, 0.012658, 0.0, 0.0, 0.012658],
    #        help="matdyn qpoints")

    # ------------------------------------------------------------- 
    # for server
    # -------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.emuhelper/server.conf")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

    matdyn_input["asr"] = args.asr

    #matdyn_qpoints = []
    #for i in range(0, len(args.matdyn_qpoints), 4):
    #    matdyn_qpoints.append([float(args.matdyn_qpoints[i]), float(args.matdyn_qpoints[i+1]), float(args.matdyn_qpoints[i+2]), float(args.matdyn_qpoints[i+3])])

    #print(matdyn_qpoints)

    task = dfpt_run()
    task.get_xyz(args.file)
    task.set_matdyn(matdyn_input=matdyn_input)
    task.matdyn(directory=args.directory, mpi=args.mpi, runopt=args.runopt) # nqpoints=args.nqpoints, qpoints=matdyn_qpoints)

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
        ctl.submit(workdir=args.directory, jobfile="matdyn.sub")
