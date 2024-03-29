#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.dfpt import dfpt_run


"""
usage:
"""

matdyn_input = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-static",
            help="directory for the static running")

    parser.add_argument("-f", "--file", type=str,
            help="the xyz file")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

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
    #                       run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="qe-matdyn",
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

    matdyn_input["asr"] = args.asr

    # get qpoints either from --qpoints or from --qpoints-file
    # qpoints is in format like this:
    # [[kx, ky, kz, xcoord, label], ...] like [[0.0, 0,0, 0.0, 0.0, 'GAMMA']]
    # if the label is a str like 'GAMMA', 'K', etc, the q point is a specialk,
    # if the label is None, then the q points is not a special
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
    task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)     
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
    task.matdyn(directory=args.directory, runopt=args.runopt, auto=args.auto)
