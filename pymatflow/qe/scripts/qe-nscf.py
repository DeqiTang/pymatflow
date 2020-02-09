#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.static import static_run
from pymatflow.remote.server import server_handle

"""
usage:
"""

control = {}
system = {}
electrons = {}
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-static",
            help="Directory for the static running.")
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

    parser.add_argument("--ecutrho", type=int, default=None,
            help="Kinetic energy cutoff for charge density and potential in unit of Rydberg, default value: None")

    parser.add_argument("--kpoints-option", type=str, default="automatic", 
            choices=["automatic", "gamma", "crystal_b"],
            help="Kpoints generation scheme option for the SCF or non-SCF calculation")

    parser.add_argument("--kpoints-mp", type=int, nargs="+",
            default=[4, 4, 4, 0, 0, 0],
            help="Monkhorst-Pack kpoint grid, in format like --kpoints-mp 4 4 4 0 0 0")

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

    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")
    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")
    parser.add_argument("--jobname", type=str, default="pwscf-nscf",
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
    electrons["conv_thr"] = args.conv_thr

    task = static_run()
    task.get_xyz(xyzfile)
    task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
    task.set_params(control=control, system=system, electrons=electrons)
    task.nscf(directory=args.directory, runopt=args.runopt, mpi=args.mpi, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)

    server_handle(auto=args.auto, directory=args.directory, jobfilebase="static-nscf", server=args.server)
