#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from pymatflow.qe.static import static_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync

"""
usage qe-converge-kpoints.py -f xxx.xyz --range nk_min nk_max step
"""

control_params = {}
system_params = {}
electrons_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of the calculation", type=str, default="tmp-qe-ecutrho")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--runopt", help="gen run genrun", type=str, default="genrun")
    parser.add_argument("--ecutwfc", help="better specify a converged ecutwfc", type=int, default=100)
    parser.add_argument("--ecutrho", help="better specify a converged ecutrho", type=int, default=400)
    parser.add_argument("--range", help="nk_min nk_max step", nargs='+', type=int)
    parser.add_argument("--conv-thr", help="conv_thr", type=float, default=1.0e-6)

    parser.add_argument("--occupations", type=str, default="smearing",
            choices=["smearing", "tetrahedra", "tetrahedra_lin", "tetrahedra_opt", "fixed", "from_input"],
            help="Occupation method for the calculation.")
    
    parser.add_argument("--smearing", type=str, default="gaussian",
            choices=["gaussian", "methfessel-paxton", "marzari-vanderbilt", "fermi-dirac"],
            help="Smearing type for occupations by smearing, default is gaussian in this script")

    parser.add_argument("--degauss", type=float, default=0.001,
            help="Value of the gaussian spreading (Ry) for brillouin-zone integration in metals.(defualt: 0.001 Ry)")

    parser.add_argument("--vdw-corr", type=str, default="none",
            choices=["dft-d", "dft-d3", "ts", "xdm"],
            help="Type of Van der Waals correction in the calculation")

    # for server
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.emuhelper/server.conf")
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file
    system_params["ecutwfc"] = args.ecutwfc
    system_params["ecutrho"] = args.ecutrho
    system_params["occupations"] = args.occupations
    system_params["smearing"] = args.smearing
    system_params["degauss"] = args.degauss
    system_params["vdw_corr"] = args.vdw_corr
    electrons_params["conv_thr"] = args.conv_thr

    task = static_run()
    task.get_xyz(xyzfile)
    task.set_params(control=control_params, system=system_params, electrons=electrons_params)
    task.converge_kpoints(args.range[0], args.range[1], args.range[2], directory=args.directory, runopt=args.runopt)

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
        ctl.submit(workdir=args.directory, jobfile="converge-kpoints.sub")
