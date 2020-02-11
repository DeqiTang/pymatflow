#!/usr/bin/env python
 # _*_ coding: utf-8 _*_

import argparse

from pymatflow.siesta.opt import opt_run
from pymatflow.remote.server import server_handle

"""
usage:
   siesta-opt.py xxx.xyz
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # ===========================
    # general parameters
    # ===========================
    parser.add_argument("-d", "--directory", type=str, default="tmp-siesta-opt",
            help="directory to do the optimization calculation")

    parser.add_argument("-f", "--file", type=str,
            help="the xyz structure file with second line specifying cell parameters")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command")

    #parser.add_argument("-p", "--properties" ,help="Option for properties calculation", type=int, default=0)
    parser.add_argument("-m", "--mode", type=int, default=0,
            choices=[0, 1],
            help="Optimization mode: 0(not-variable-cell), 1(variable-cell)")

    # =========================
    #      Electrons
    # =========================
    parser.add_argument("--meshcutoff", type=int, default=200,
            help="MeshCutoff (Ry)")

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
            help="OccupationFunction(FD or MP)")

    parser.add_argument("--electronic-temperature", type=int ,default=300,
            help="Electronic Temperature")

    # ==================================================
    #           ions relaed parameter
    # ==================================================
    parser.add_argument("--vc", type=str, default="false",
            choices=["true", "false"],
            help="MD.VariableCell")

    parser.add_argument("--forcetol", type=float, default=0.04,
            help="Force tolerance in coordinate optimization. default=0.04 eV/Ang")

    parser.add_argument("--stresstol", type=float, default=1,
            help="Stress tolerance in variable-cell CG optimization. default=1 GPa")

    parser.add_argument("--targetpressure", type=float, default=0,
            help="Target pressure for Parrinello-Rahman method, variable cell optimizations, and annealing options.")

    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="siesta-opt",
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
    directory = args.directory

    params = {}

    params["MeshCutoff"] = args.meshcutoff
    params["SolutionMethod"] = args.solution_method
    params["XC.funtional"] = args.functional
    params["XC.authors"] = args.authors
    params["DM.Tolerance"] = args.tolerance
    params["DM.NumberPulay"] = args.numberpulay
    params["DM.MixingWeight"] = args.mixing
    params["OccupationFunction"] = args.occupation
    params["ElectronicTemperature"] = args.electronic_temperature

    params["MD.VariableCell"] = args.vc
    params["MD.MaxForceTol"] = args.forcetol
    params["MD.MaxStressTol"] = args.stresstol
    params["MD.TargetPressure"] = args.targetpressure

    task = opt_run()
    task.get_xyz(args.file)
    task.set_params(params=params)
    task.set_kpoints(kpoints_mp=args.kpoints_mp)
    task.opt(directory=directory, runopt=args.runopt, mpi=args.mpi, mode=args.mode)

    # server handle
    server_handle(auto=args.auto, directory=args.directory, jobfilebase="geometric-optimization", server=args.server)
