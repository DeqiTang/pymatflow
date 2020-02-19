#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.static import static_run


"""
usage:
    qe-bands.py -f xxx.xyz
"""

control = {}
system = {}
electrons = {}
bands = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-static",
            help="Directory for the static running.")

    parser.add_argument("-f", "--file", type=str,
            help="The xyz file containing the structure to be simulated.")

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
            help="Kinetic energy cutoff for charge density and potential in unit of Rydberg, default value: None")

    parser.add_argument("--kpoints-option", type=str, default="crystal_b",
            choices=["automatic", "gamma", "crystal_b"],
            help="Kpoints generation scheme option for band calculation")

    parser.add_argument("--kpoints-mp", type=str, nargs="+",
            default=[1, 1, 1, 0, 0, 0],
            help="Monkhorst-Pack kpoint grid, in format like --kpoints-mp 1 1 1 0 0 0")

    parser.add_argument("--crystal-b", type=str, nargs="+", default=None,
            help="manual input kpath in crystal_b, like --crystal-b '0.000000 0.000000 0.000000 GAMMA 5' '0.500000 0.000000 0.000000 X 5' '0.0000 0.000 0.50000 A |' '0.5 0.5 0.5 R '")

    parser.add_argument("--crystal-b-file", type=str, default='kpath-from-seekpath.txt',
            help="manual input kpath in crystal_b read from the file")

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

    # -----------------------------------------
    #         bands.x related parameters
    # -----------------------------------------
    parser.add_argument("--lsym", type=str, default=".true.",
            choices=[".true.", ".false."],
            help="set lsym variable in bands.x input.")

    # -----------------------------------------------------------------
    #                       run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="qe-band-structure",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

    system["occupations"] = args.occupations
    system["smearing"] = args.smearing
    system["degauss"] = args.degauss

    bands["lsym"] = args.lsym

    # --------------------------------------------------------------
    # process crystal_b

    if args.crystal_b != None:
        # crystal_b from script argument args.crystal_b
        crystal_b = []
        for kpoint in args.crystal_b:
            if kpoint.split()[4] != "|":
                crystal_b.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    int(kpoint.split()[4]),
                    ])
            elif kpoint.split()[4] == "|":
                crystal_b.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    "|",
                    ])
    elif args.crystal_b == None:
        # crystal_b read from file specified by args.crystal_b_file
        # file is in format like this
        """
        5
        0.0 0.0 0.0 #GAMMA 15
        x.x x.x x.x #XXX |
        x.x x.x x.x #XXX 10
        x.x x.x x.x #XXX 15
        x.x x.x x.x #XXX 20
        """
        # if there is a '|' behind the label it means the path is
        # broken after that point!!!
        crystal_b = []
        with open(args.crystal_b_file, 'r') as fin:
            crystal_b_file = fin.readlines()
        nk = int(crystal_b_file[0])
        for i in range(nk):
            if crystal_b_file[i+1].split("\n")[0].split()[4] != "|":
                crystal_b.append([
                    float(crystal_b_file[i+1].split()[0]),
                    float(crystal_b_file[i+1].split()[1]),
                    float(crystal_b_file[i+1].split()[2]),
                    crystal_b_file[i+1].split()[3].split("#")[1].upper(),
                    int(crystal_b_file[i+1].split()[4]),
                    ])
            elif crystal_b_file[i+1].split("\n")[0].split()[4] == "|":
                crystal_b.append([
                    float(crystal_b_file[i+1].split()[0]),
                    float(crystal_b_file[i+1].split()[1]),
                    float(crystal_b_file[i+1].split()[2]),
                    crystal_b_file[i+1].split()[3].split("#")[1].upper(),
                    '|',
                    ])
    else:
        pass
    # --------------------------------------------------------------------


    task = static_run()
    task.get_xyz(args.file)
    task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp, crystal_b=crystal_b)
    task.set_params(control=control, system=system, electrons=electrons)
    task.set_bands(bands_input=bands)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.bands(directory=args.directory, runopt=args.runopt, auto=args.auto)
