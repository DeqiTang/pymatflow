#!/usr/bin/env python
 # _*_ coding: utf-8 _*_

import argparse

from emuhelper.siesta.opt import opt_run

"""
usage:
   siesta-opt.py xxx.xyz
"""
electrons = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # ===========================
    # general parameters
    # ===========================
    parser.add_argument("-d", "--directory", type=str, default="tmp-siesta-opt"
            help="directory of the calculation")
    parser.add_argument("-f", "--file", type=str,
            help="the xyz file name")
    parser.add_argument("--mpi", type=str, default="",
            help="MPI command")
    #parser.add_argument("-p", "--properties" ,help="Option for properties calculation", type=int, default=0)
    parser.add_argument("-m", "--mode", type=int, default=0,
            choices=[0, 1]
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
    parser.add_argument("-k", "--kpoints", type=str, default="3 3 3",
            help="set kpoints like '3 3 3'")
    parser.add_argument("--occupation", type=str, default="FD",
            help="OccupationFunction(FD or MP)")
    parser.add_argument("--electronic-temperature", type=int ,default=300,
            help="Electronic Temperature")
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()
    xyzfile = args.file
    directory = args.directory
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(3)]
    
    electrons["MeshCutoff"] = args.meshcutoff
    electrons["SolutionMethod"] = args.solution_method
    electrons["XC.funtional"] = args.functional
    electrons["XC.authors"] = args.authors
    electrons["DM.Tolerance"] = args.tolerance
    electrons["DM.NumberPulay"] = args.numberpulay
    electrons["DM.MixingWeight"] = args.mixing
    electrons["OccupationFunction"] = args.occupation
    electrons["ElectronicTemperature"] = args.electronic_temperature

    task = opt_run(xyzfile)
    task.opt(directory=directory, runopt="genrun", mpi=args.mpi, electrons=electrons, kpoints_mp=kpoints_mp, mode=args.mode)
    task.analysis()
