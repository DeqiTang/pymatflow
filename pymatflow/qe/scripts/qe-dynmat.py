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

dynmat_input = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()    
    parser.add_argument("-d", "--directory", help="directory for the static running", type=str, default="tmp-qe-static")
    parser.add_argument("--mpi", help="MPI commadn", type=str, default="")
    parser.add_argument("-f", "--file", help="the xyz file", type=str)
    parser.add_argument("--runopt", help="gen, run, or genrun", type=str, default="genrun")
   
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
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")
    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"]
            help="type of remote server, can be pbs or yh")

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
    task.dynmat(directory=args.directory, mpi=args.mpi, runopt=args.runopt)

    # server handle
    if args.auto == 0:
        pass
    elif args.auto == 1:
        mover = rsync()
        if args.server == "pbs":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_pbs.conf"))
            pass
        elif args.server == "yh":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_yh.conf"))
        mover.copy_default(source=os.path.abspath(args.directory))
    elif args.auto == 2:
        mover = rsync()
        if args.server == "pbs":
            pass
        elif args.server == "yh":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_yh.conf"))
        mover.copy_default(source=os.path.abspath(args.directory))
        ctl = ssh()
        if args.server == "pbs":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_pbs.conf"))
            ctl.login()
            ctl.submit(workdir=args.directory, jobfile="dynmat.pbs", server="pbs")
        elif args.server == "yh":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_yh.conf"))
            ctl.login()
            ctl.submit(workdir=args.directory, jobfile="dynmat.sub", server="yh")
