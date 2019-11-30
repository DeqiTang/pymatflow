#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.opt import opt_run
from pymatflow.remote.rsync import rsync
from pymatflow.remote.ssh import ssh
"""
usage: qe-vc-relax.py -f xxx.xyz
"""

control_params = {}
system_params = {}
electrons_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", 
            type=str, default="tmp-qe-vc-relax",
            help="directory for the vc-relax running")

    parser.add_argument("-f", "--file", 
            type=str,
            help="the xyz file containg the structure to be simulated")

    parser.add_argument("--runopt",
            type=str, default="genrun",
            help="run option, could be: gen, run, genrun")

    parser.add_argument("--mpi", 
            type=str, default="",
            help="the mpi command used")

    # -------------------------------------------------------------------
    #                       scf related parameters
    # -------------------------------------------------------------------
    parser.add_argument("--ecutwfc", 
            type=int, default=100)

    parser.add_argument("--ecutrho", 
            type=int, default=400)

    parser.add_argument("--kpoints-option", help="kpoints option", type=str, default="automatic")

    parser.add_argument("-k", "--kpoints", help="set kpoints like '1 1 1 0 0 0'", type=str, default="1 1 1 0 0 0")

    parser.add_argument("--conv-thr", type=float, default=1.0e-6,
            help="the conv_thr for scf, when doing geometric optimization better use a strict covnergec for scf")

    parser.add_argument("--occupations", help="occupation type", type=str, default="smearing")

    parser.add_argument("--smearing", help="smearing type", type=str, default="gaussian")

    parser.add_argument("--degauss", help="value of the gaussian spreading (Ry) for brillouin-zone integration in metals.", type=float, default=0.001)
    parser.add_argument("--vdw-corr", help="vdw_corr = dft-d, dft-d3, ts-vdw, xdm", type=str, default="none")


    # -------------------------------------------------------------------
    #               geometric optimization related parameters
    # -------------------------------------------------------------------
    parser.add_argument("--etot-conv-thr", 
            type=float, default=1.0e-4,
            help="convergence threshold of energy for geometric optimization")

    parser.add_argument("--forc-conv-thr", 
            type=float, default=1.0e-3,
            help="convergence threshold for force in optimization,(usually it is more important than energy)")

    parser.add_argument("--nstep",
            type=int, default=50,
            help="maximum ion steps for geometric optimization")
    # -------------------------------------------------------------------
    #                       for server handling
    # -------------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.emuhelper/server.conf")


    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file
    control_params["etot_conv_thr"] = args.etot_conv_thr
    control_params["forc_conv_thr"] = args.forc_conv_thr
    control_params["nstep"] = args.nstep
    system_params["ecutwfc"] = args.ecutwfc
    system_params["ecutrho"] = args.ecutrho
    system_params["occupations"] = args.occupations
    system_params["smearing"] = args.smearing
    system_params["degauss"] = args.degauss
    system_params["vdw_corr"] = args.vdw_corr
    electrons_params["conv_thr"] = args.conv_thr
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]

    task = opt_run(xyzfile)
    task.vc_relax(directory=args.directory, runopt=args.runopt, mpi=args.mpi, control=control_params, system=system_params, electrons=electrons_params, kpoints_option=args.kpoints_option, kpoints_mp=kpoints_mp)

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
        ctl.submit(workdir=args.directory, jobfile="vc-relax.in.sub")
