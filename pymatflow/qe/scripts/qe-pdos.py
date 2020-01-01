#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.static import static_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync

"""
usage:
    qe-pdos.py -f xxx.xyz
"""

projwfc_input = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()    
    parser.add_argument("-d", "--directory", help="directory for the static running", type=str, default="tmp-qe-static")
    parser.add_argument("--mpi", help="MPI commadn", type=str, default="")
    parser.add_argument("-f", "--file", help="the xyz file", type=str)
    parser.add_argument("--runopt", help="gen, run, or genrun", type=str, default="genrun")
    parser.add_argument("--filpdos", help="output projected dos file name", type=str, default="projwfc")
    parser.add_argument("--ngauss", help="gaussian broadening type", type=str, default='default')
    parser.add_argument("--degauss", help="gaussian broadening", type=str, default='default')
    parser.add_argument("--emin", help="min energy for DOS", type=str, default='default')
    parser.add_argument("--emax", help="max energy for DOS", type=str, default='default')
    parser.add_argument("--deltae", help="DeltaE: energy grid step (eV)", type=str, default='default')

    # for server
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.emuhelper/server.conf")
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

    if args.ngauss == 'default':
        ngauss = args.ngauss
    else:
        ngauss = int(args.ngauss)
    if args.degauss == 'default':
        degauss = args.degauss
    else:
        degauss = float(args.degauss)
    if args.emin == 'default':
        emin = args.emin
    else:
        emin = float(args.emin)
    if args.emax == 'default':
        emax = args.emax
    else:
        emax = float(args.emax)
    if args.deltae == 'default':
        deltae = args.deltae
    else:
        deltae = float(args.deltae)

    projwfc_input["filpdos"] = args.filpdos
    projwfc_input["ngauss"] = ngauss
    projwfc_input["degauss"] = degauss
    projwfc_input["emin"] = emin
    projwfc_input["emax"] = emax
    projwfc_input["deltae"] = deltae

    task = static_run()
    task.get_xyz(args.file)
    task.set_projwfc(projwfc_input=projwfc_input)
    task.projwfc(directory=args.directory, runopt=args.runopt, mpi=args.mpi)

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
        ctl.submit(workdir=args.directory, jobfile="static-projwfc.sub")
