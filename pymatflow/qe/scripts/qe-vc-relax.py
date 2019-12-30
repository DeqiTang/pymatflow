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

control = {}
system = {}
electrons = {}
ions = {}
cell = {}

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

    parser.add_argument("--kpoints-option", type=str, default="automatic", 
            choices=["automatic", "gamma", "tpiba_b"],
            help="Kpoints generation scheme option for the SCF or non-SCF calculation")

    parser.add_argument("--kpoints-mp", type=int, nargs="+",
            default=[1, 1, 1, 0, 0, 0],
            help="Monkhorst-Pack kpoint grid, in format like --kpoints-mp 1 1 1 0 0 0")

    parser.add_argument("--conv-thr", type=float, default=1.0e-6,
            help="the conv_thr for scf, when doing geometric optimization better use a strict covnergec for scf")

    parser.add_argument("--occupations", type=str, default="smearing",
            choices=["smearing", "tetrahedra", "tetrahedra_lin", "tetrahedra_opt", "fixed", "from_input"],
            help="Occupation method for the calculation.")
    
    parser.add_argument("--smearing", type=str, default="gaussian",
            choices=["gaussian", "methfessel-paxton", "marzari-vanderbilt", "fermi-dirac"],
            help="Smearing type for occupations by smearing, default is gaussian in this script")

    parser.add_argument("--degauss", type=float, default=0.001,
            help="Value of the gaussian spreading (Ry) for brillouin-zone integration in metals.(defualt: 0.001 Ry)")

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
    control["etot_conv_thr"] = args.etot_conv_thr
    control["forc_conv_thr"] = args.forc_conv_thr
    control["nstep"] = args.nstep
    system["ecutwfc"] = args.ecutwfc
    system["ecutrho"] = args.ecutrho
    system["occupations"] = args.occupations
    system["smearing"] = args.smearing
    system["degauss"] = args.degauss
    system["vdw_corr"] = args.vdw_corr
    electrons["conv_thr"] = args.conv_thr

    task = opt_run()
    task.get_xyz(xyzfile)
    task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
    task.set_params(control=control, system=system, electrons=electrons, ions=ions, cell=cell)
    task.vc_relax(directory=args.directory, runopt=args.runopt, mpi=args.mpi)

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
