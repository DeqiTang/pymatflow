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
    parser.add_argument("-d", "--directory", help="directory of the calculation", type=str, default="tmp-siesta-opt")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--mpi", help="MPI command", type=str, default="")
    #parser.add_argument("-p", "--properties" ,help="Option for properties calculation", type=int, default=0)
    parser.add_argument("-m", "--mode", help="Optimization mode: 0(no vc), 1(vc)", type=int, default=0)
    parser.add_argument("--meshcutoff", help="MeshCutoff (Ry)", type=int, default=200)
    parser.add_argument("--solution-method", help="SolutionMethod(diagon, OMM, OrderN, PEXSI)", type=str, default="diagon")
    parser.add_argument("--functional", help="XC.functional", type=str, default="GGA")
    parser.add_argument("--authors", help="XC.authors", type=str, default="PBE")
    parser.add_argument("--tolerance", help="DM.Tolerance", type=float, default=1.0e-6)
    parser.add_argument("--numberpulay", help="DM.NumberPulay", type=int ,default=8)
    parser.add_argument("--mixing", help="DM.MixingWeight", type=float, default=0.1)
    parser.add_argument("-k", "--kpoints", help="set kpoints like '3 3 3'", type=str, default="3 3 3")
    parser.add_argument("--occupation", help="OccupationFunction(FD or MP)", type=str, default="FD")
    parser.add_argument("--electronic-temperature", help="Electronic Temperature", type=int, default=300)
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
