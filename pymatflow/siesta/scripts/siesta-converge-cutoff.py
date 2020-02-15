#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from pymatflow.siesta.static import static_run
from pymatflow.remote.server import server_handle

"""
usage:
    siesta-converge-cutoff.py xxx.xyz emin emax step
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-siesta-converge-cutoff",
            help="directory to do the MeshCutoff converge test calculation")

    parser.add_argument("-f", "--file", type=str,
            help="the xyz structure file with second line specifying cell parameters")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--mpi", help="MPI command", default="")

    parser.add_argument("--range", type=int, nargs="+",
            help="Test Range of MeshCutoff (Ry)")

    parser.add_argument("--solution-method", type=str, default="diagon",
            help="SolutionMethod(diagon, OMM, OrderN, PEXSI)")

    parser.add_argument("--functional", type=str, default="GGA",
            help="XC.functional")

    parser.add_argument("--authors", type=str, default="PBE",
            help="XC.authors")

    parser.add_argument("--tolerance", type=float, default=1.0e-6,
            help="DM.Tolerance")

    parser.add_argument("--numberpulay", type=int, default=8,
            help="DM.NumberPulay")

    parser.add_argument("--mixing", type=float, default=0.1,
            help="DM.MixingWeight")

    parser.add_argument("--kpoints-mp", type=int, nargs="+",
            default=[3, 3, 3],
            help="set kpoints like '3 3 3'")

    parser.add_argument("--occupation", type=str, default="FD",
            chioces=["FD", "MP"],
            help="OccupationFunction(FD or MP)",)

    parser.add_argument("--electronic-temperature", type=int, default=300,
            help="Electronic Temperature")


    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="siesta-scf",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file

    params = {}

    params["SolutionMethod"] = args.solution_method
    params["XC.funtional"] = args.functional
    params["XC.authors"] = args.authors
    params["DM.Tolerance"] = args.tolerance
    params["DM.NumberPulay"] = args.numberpulay
    params["DM.MixingWeight"] = args.mixing
    params["OccupationFunction"] = args.occupation
    params["ElectronicTemperature"] = args.electronic_temperature

    task = static_run()
    task.get_xyz(args.file)
    task.set_params(params=params)
    task.set_kpoints(kpoints_mp=args.kpoints_mp)
    task.converge_cutoff(args.range[0], args.range[1], args.range[2], directory=args.directory, runopt=args.runopt, mpi=args.mpi)

    # server handle
    server_handle(auto=args.auto, directory=args.directory, jobfilebase="converge-cutoff", server=args.server)
