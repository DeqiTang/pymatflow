#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import argparse

from pymatflow.qe.md import md_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync

"""
usage: qe-md.py xxx.xyz
"""

control = {}
system = {}
electrons = {}
ions = {}
cell = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory for the vc-md running", type=str, default="tmp-qe-vc-md")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--runopt", help="gen, run, or genrun", type=str, default="genrun")
    parser.add_argument("--mpi", help="MPI command", type=str, default="")
    parser.add_argument("--nstep", help="maximum ion steps", type=int, default=50)
    parser.add_argument("--ecutwfc", help="ecutwfc", type=int, default=100)
    parser.add_argument("--ecutrho", help="ecutrho", type=int, default=400)
    
    parser.add_argument("--kpoints-option", type=str, default="automatic", 
            choices=["automatic", "gamma", "tpiba_b"],
            help="Kpoints generation scheme option for the SCF or non-SCF calculation")

    parser.add_argument("--kpoints-mp", type=int, nargs="+",
            default=[1, 1, 1, 0, 0, 0],
            help="Monkhorst-Pack kpoint grid, in format like --kpoints-mp 1 1 1 0 0 0")
    
    parser.add_argument("--conv-thr", help="conv_thr of scf", type=float, default=1.e-6)

    parser.add_argument("--occupations", type=str, default="smearing",
            choices=["smearing", "tetrahedra", "tetrahedra_lin", "tetrahedra_opt", "fixed", "from_input"],
            help="Occupation method for the calculation.")
    
    parser.add_argument("--smearing", type=str, default="gaussian",
            choices=["gaussian", "methfessel-paxton", "marzari-vanderbilt", "fermi-dirac"],
            help="Smearing type for occupations by smearing, default is gaussian in this script")

    parser.add_argument("--degauss", type=float, default=0.001,
            help="Value of the gaussian spreading (Ry) for brillouin-zone integration in metals.(defualt: 0.001 Ry)")

    parser.add_argument("--vdw-corr", help="vdw_corr = dft-d, dft-d3, ts-vdw, xdm", type=str, default='none')

    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")
    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"]
            help="type of remote server, can be pbs or yh")

    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
    # ==========================================================   

    args = parser.parse_args()
    xyzfile = args.file
    control["nstep"] = args.nstep
    system["ecutwfc"] = args.ecutwfc
    system["ecutrho"] = args.ecutrho
    system["occupations"] = args.occupations
    system["smearing"] = args.smearing
    system["degauss"] = args.degauss
    system["vdw_corr"] = args.vdw_corr
    electrons["conv_thr"] = args.conv_thr
 
    task = md_run()
    task.get_xyz(args.file)
    task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
    task.set_params(control=control, system=system, electrons=electrons, ions=ions, cell=cell)
    task.vc_md(directory=args.directory, runopt=args.runopt, mpi=args.mpi)

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
            ctl.submit(workdir=args.directory, jobfile="vc-md.pbs", server="pbs")
        elif args.server == "yh":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_yh.conf"))
            ctl.login()
            ctl.submit(workdir=args.directory, jobfile="vc-md.sub", server="yh")
