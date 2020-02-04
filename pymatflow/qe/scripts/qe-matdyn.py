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
    parser.add_argument("--qpoints", type=str, nargs="+", default=None,
            help="matdyn qpoints manual input like --qpoints '0.0 0.0 0.0 0.0 GAMMA' '0.5 0.0 0.0 0.5' 'xxx' 'xxx', if the qpoints is not a special q point, do not specify the label!!!!!!!!!!!!!!!!!!")
    parser.add_argument("--qpoints-file", type=str, default="matdyn-qpoints.txt",
            help="file to get the qpoints for matdyn")

    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")
    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")
    parser.add_argument("--jobname", type=str, default="qe-matdyn",
            help="jobname on the pbs server")
    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")
    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

    matdyn_input["asr"] = args.asr

    # get qpoints either from --qpoints or from --qpoints-file
    qpoints = []
    if args.qpoints != None:
       # qpoints read from script argument args.qpoints
       for point in args.qpoints:
           if len(point.split()) == 4:
               qpoints.append([
                   float(point.split()[0]), 
                   float(point.split()[1]), 
                   float(point.split()[2]),
                   float(point.split()[3]),
                   None,
                   ])
           elif len(point.split()) == 5:
                qpoints.append([
                   float(point.split()[0]), 
                   float(point.split()[1]), 
                   float(point.split()[2]),
                   float(point.split()[3]),
                   point.split()[4].upper(),
                   ])
    elif args.qpoints == None:
        # qpoints read from a file contains manual qpoints setting
        # qpoints_file in format like this:
        """
        3
        0.000000 0.000000 0.000000 0.000000 #GAMMA
        0.500000 0.000000 0.000000 0.500000
        0.000000 0.500000 0.000000 0.500000
        """
        # if a q point in qpoints_file is a special q there must be a
        # comment after that line to give the special q label
        #
        with open(args.qpoints_file, 'r') as fin:
          lines = fin.readlines()
        nq = int(lines[0])
        for i in range(nq):
            if len(lines[i+1].split()) == 0:
                continue
            if len(lines[i+1].split()) == 4:
                qpoints.append([
                    float(lines[i+1].split()[0]), 
                    float(lines[i+1].split()[1]), 
                    float(lines[i+1].split()[2]),
                    float(lines[i+1].split()[3]),
                    None,
                ])
            else:
                qpoints.append([
                    float(lines[i+1].split()[0]), 
                    float(lines[i+1].split()[1]), 
                    float(lines[i+1].split()[2]),
                    float(lines[i+1].split()[3]),
                    lines[i+1].split("\n")[0].split("#")[1].upper(),
                ])
        #
        # get information on high symmetry q point



    task = dfpt_run()
    task.get_xyz(args.file)
    task.set_matdyn(matdyn_input=matdyn_input, qpoints=qpoints)
    task.matdyn(directory=args.directory, mpi=args.mpi, runopt=args.runopt, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)

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
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_pbs.conf"))
        elif args.server == "yh":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_yh.conf"))
        mover.copy_default(source=os.path.abspath(args.directory))
        ctl = ssh()
        if args.server == "pbs":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_pbs.conf"))
            ctl.login()
            ctl.submit(workdir=args.directory, jobfile="matdyn.pbs", server="pbs")
        elif args.server == "yh":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_yh.conf"))
            ctl.login()
            ctl.submit(workdir=args.directory, jobfile="matdyn.sub", server="yh")
