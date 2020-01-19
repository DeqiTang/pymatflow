#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.neb import neb_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync

"""
usage:
    qe-neb.py -f xxx.xyz -k '2 2 2 0 0 0' --ecutwfc 100
"""


control = {}
system = {}
electrons = {}
path = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of the calculation", type=str, default="tmp-qe-neb")
    parser.add_argument("--restart-mode", help="restart_mode", type=str, default="from_scratch")
    parser.add_argument("--runopt", help="gen, run, or genrun", type=str, default="genrun")
    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")
    parser.add_argument("--images", help="the image xyz file(--images first.xyz imtermediate-1.xyz intermediate-2.xyz ... last.xyz)", nargs='+', type=str)


    # -------------------------------------------------------------------
    #                       scf related parameters
    # -------------------------------------------------------------------
    parser.add_argument("--ecutwfc", type=int, default=100,
            help="Kinetic energy cutoff for wave functions in unit of Rydberg, default value: 100 Ry")

    parser.add_argument("--ecutrho", type=int, default=400,
            help="Kinetic energy cutoff for charge density and potential in unit of Rydberg, default value: 400 Ry")

    parser.add_argument("--kpoints-option", type=str, default="automatic", 
            choices=["automatic", "gamma", "tpiba_b"],
            help="Kpoints generation scheme option for the SCF or non-SCF calculation")

    parser.add_argument("--kpoints-mp", type=int, nargs="+",
            default=[1, 1, 1, 0, 0, 0],
            help="Monkhorst-Pack kpoint grid, in format like --kpoints-mp 1 1 1 0 0 0")

    parser.add_argument("--conv-thr", type=float, default=1.0e-6,
            help="Convergence threshold for SCF calculation.")

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

    parser.add_argument("--nbnd", type=int, default=None,
            help="Number of electronic states (bands) to be calculated")

    # --------------------------------------------------------------
    #                params for neb namelist &path
    # --------------------------------------------------------------
    parser.add_argument("--string-method", help="string_method", type=str, default="neb")
    parser.add_argument("--nstep-path", help="nstep_path", type=int, default=100)
    parser.add_argument("--opt-scheme", help="Specify the type of optimization scheme(sd, broyden, broyden2, quick-min, langevin)", type=str, default="broyden")

    parser.add_argument("--num-of-images", help="number of total images(including the initial and final image). about how to set proper number of images: usually the inter-image distance between 1~2Bohr is OK", type=int, default=5)

    parser.add_argument("--k-max", help="Set them to use a Variable Elastic Constants scheme elastic constants are in the range [ k_min, k_max  ], this is useful to rise the resolution around the saddle point", type=float, default=0.3e0)
    parser.add_argument("--k-min", help="Set them to use a Variable Elastic Constants scheme elastic constants are in the range [ k_min, k_max  ], this is useful to rise the resolution around the saddle point", type=float, default=0.2e0)
    parser.add_argument("--ci-scheme", help="Specify the type of Climbing Image scheme(no-CI, auto, manual)", type=str, default="auto")
    parser.add_argument("--path_thr", help="path_thr", type=float, default=0.05)
    parser.add_argument("--ds", help="Optimisation step length ( Hartree atomic units )", type=float, default=1.e0)
    parser.add_argument("--first-last-opt", type=bool, default=False)
    
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
    directory = args.directory
    
    system["ecutwfc"] = args.ecutwfc
    system["occupations"] = args.occupations
    system["smearing"] = args.smearing
    system["degauss"] = args.degauss
    system["vdw_corr"] = args.vdw_corr
    electrons["conv_thr"] = args.conv_thr

    path["string_method"] = args.string_method
    path["nstep_path"] = args.nstep_path
    path["opt_scheme"] = args.opt_scheme
    path["num_of_images"] = args.num_of_images
    path["k_max"] = args.k_max
    path["k_min"] = args.k_min
    path["CI_scheme"] = args.ci_scheme
    path["path_thr"] = args.path_thr
    path["ds"] = args.ds
    path["first_last_opt"] = args.first_last_opt

    task = neb_run()
    task.get_images(images=args.images)
    task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
    task.set_path(path=path)
    task.neb(directory=directory, mpi=args.mpi, runopt=args.runopt, restart_mode=args.restart_mode)

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
            ctl.submit(workdir=args.directory, jobfile="neb.pbs", server="pbs")
        elif args.server == "yh":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_yh.conf"))
            ctl.login()
            ctl.submit(workdir=args.directory, jobfile="neb.sub", server="yh")
