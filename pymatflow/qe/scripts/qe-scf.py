#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.static import static_run

"""
usage:
    qe-scf.py -f xxx.xyz -k '2 2 2 0 0 0' --ecutwfc 100
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

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")


    # -------------------------------------------------------------------
    #                       scf related parameters
    # -------------------------------------------------------------------
    parser.add_argument("--ecutwfc", type=int, default=100,
            help="Kinetic energy cutoff for wave functions in unit of Rydberg, default value: 100 Ry")

    parser.add_argument("--ecutrho", type=int, default=None,
            help="Kinetic energy cutoff for charge density and potential in unit of Rydberg, default vaslue: None")

    parser.add_argument("--kpoints-option", type=str, default="automatic",
            choices=["automatic", "gamma", "crystal_b"],
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

    parser.add_argument("--tstress", type=str, default=".false.",
            choices=[".true.", ".false."],
            help="calculate stress. default=.false.")

    # magnetic related parameters
    parser.add_argument("--nspin", type=int, default=None,
            choices=[1, 2],
            help="choose either 1 or 2, and 4 should not be used as suggested by pwscf official documentation.")

    parser.add_argument("--starting-magnetization", type=float, nargs="+", default=None,
            help="starting_magnetization(i), i=1,ntyp -> Starting spin polarization on atomic type i in a spin polarized calculation. Values range between -1 (all spins down for the valence electrons of atom type i) to 1 (all spins up).")

    parser.add_argument("--noncolin", type=str, default=None,
            choices=[".true.", ".false."],
            help="if .true. the program will perform a noncollinear calculation.")
    # -----------------------------------------------------------
    #           ATOMIC_FORCES
    # -----------------------------------------------------------
    parser.add_argument("--pressure", type=float, default=None,
            help="specify pressure acting on system in unit of Pa")
    parser.add_argument("--pressuredir", type=str, default=None,
            choices=["x", "y", "z"],
            help="specify direction of pressure acting on system.")

    # -----------------------------------------------------------------
    #                       run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="pwscf-scf",
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
    control["tstress"] = args.tstress
    system["ecutwfc"] = args.ecutwfc
    system["ecutrho"] = args.ecutrho
    system["occupations"] = args.occupations
    system["smearing"] = args.smearing
    system["degauss"] = args.degauss
    system["vdw_corr"] = args.vdw_corr
    system["nbnd"] = args.nbnd
    electrons["conv_thr"] = args.conv_thr

    system["nspin"] = args.nspin
    system["starting_magnetization"] = args.starting_magnetization
    system["noncolin"] = args.noncolin


    task = static_run()
    task.get_xyz(xyzfile)
    task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
    task.set_params(control=control, system=system, electrons=electrons)
    task.set_atomic_forces(pressure=args.pressure, pressuredir=args.pressuredir)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto)
