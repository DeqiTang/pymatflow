#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.phonopy import phonopy_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync


"""
usage:
    qe-scf.py -f xxx.xyz -k '2 2 2 0 0 0' --ecutwfc 100
"""


control = {}
system = {}
electrons = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-phonopy",
            help="Directory for the phonopy running.")
    parser.add_argument("-f", "--file", type=str,
            help="The xyz file name.")
    parser.add_argument("--runopt", type=str, default="genrun", 
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")
    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    # -------------------------------------------------------------------
    #                       scf related parameters
    # -------------------------------------------------------------------
    parser.add_argument("--ecutwfc", type=int, default=100,
            help="Kinetic energy cutoff for wave functions in unit of Rydberg, default value: 100 Ry")

    parser.add_argument("--ecutrho", type=int, default=400,
            help="Kinetic energy cutoff for charge density and potential in unit of Rydberg, default value: 400 Ry")

    parser.add_argument("--kpoints-option", type=str, default="automatic", 
            choices=["automatic", "gamma", "crystal_b"],
            help="Kpoints generation scheme option for the SCF or non-SCF calculation")

    parser.add_argument("-k", "--kpoints-mp", type=int, nargs="+",
            default=[1, 1, 1, 0, 0, 0],
            help="Monkhorst-Pack kpoint grid, in format like -k 1 1 1 0 0 0")

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

    # ---------------------------------------------------------
    # Phonopy
    # ---------------------------------------------------------
    parser.add_argument("--supercell-n", type=int, nargs="+",
            default=[1, 1, 1],
            help="supercell build for Phonopy.")
    
    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")
    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")
    parser.add_argument("--jobname", type=str, default="pwscf-phonopy",
            help="jobname on the pbs server")
    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")
    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
    # ==========================================================   
    args = parser.parse_args()
    xyzfile = args.file
    system["ecutwfc"] = args.ecutwfc
    system["ecutrho"] = args.ecutrho
    system["occupations"] = args.occupations
    system["smearing"] = args.smearing
    system["degauss"] = args.degauss
    system["vdw_corr"] = args.vdw_corr
    system["nbnd"] = args.nbnd
    electrons["conv_thr"] = args.conv_thr
    

    task = phonopy_run()
    task.get_xyz(xyzfile)
    task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
    task.set_params(control=control, system=system, electrons=electrons)
    task.supercell_n = args.supercell_n
    task.phonopy(directory=args.directory, runopt=args.runopt, mpi=args.mpi, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)

    # server handle
    if args.auto == 0:
        pass
    elif args.auto == 1:
        mover = rsync()
        if args.server == "pbs":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_pbs.conf"))
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
            ctl.submit(workdir=args.directory, jobfile="phonopy-job.pbs", server="pbs")
        elif args.server == "yh":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_yh.conf"))
            ctl.login()
            ctl.submit(workdir=args.directory, jobfile="phonopy-job.sub", server="yh")
