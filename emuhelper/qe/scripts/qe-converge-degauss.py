#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.static import static_run
from emuhelper.remote.ssh import ssh
from emuhelper.remote.rsync import rsync

"""
usage qe-converge-degauss.py -f xxx.xyz --range degauss_min degauss_max step
"""

control_params = {}
# do not set occupation, smearing and degauss through system_params
system_params = {}
electrons_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--range", help="degauss_min degauss_max step", nargs='+', type=float)
    parser.add_argument("-k", "--kpoints", help="set kpoints like '1 1 1 0 0 0'", type=str, default="1 1 1 0 0 0")
    parser.add_argument("--ecutwfc", help="better a previously converged ecutwfc", type=int, default=100)

    # for server
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.emuhelper/server.conf")
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]
    system_params["ecutwfc"] = args.ecutwfc

    task = static_run(xyzfile)

    task.converge_degauss(round(args.range[0], 6), round(args.range[1], 6), round(args.range[2], 6), control=control_params, system=system_params, electrons=electrons_params, runopt="genrun")

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
        ctl.submit(workdir=args.directory, jobfile="relax.in.sub")
