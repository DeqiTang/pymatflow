#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from pymatflow.siesta.md import md_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync

"""
usage:
"""

params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # ===========================
    # general parameters
    # ===========================
    parser.add_argument("-d", "--directory", type=str, default="tmp-siesta-md",
            help="directory of the calculation")
    parser.add_argument("-f", "--file", type=str,
            help="the xyz file name")

    parser.add_argument("--runopt", type=str, default="genrun", 
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
    parser.add_argument("--mdstep", type=int, default=1000,
            help="Final time step of the MD simulation.")
    parser.add_argument("--timestep", type=float, default=1.0,
            help="Length of the time step of the MD simulation.")
    parser.add_argument("--initial-temp", type=float, default=0,
            help="Initial temperature for the MD run.")
    parser.add_argument("--target-temp", type=float, default=0,
            help="arget temperature for Nose thermostat and annealing options.")
    parser.add_argument("--vc", type=str, default="false",
            choices=["true", "false"],
            help="MD.VariableCell")


    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")
    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")
    parser.add_argument("--jobname", type=str, default="siesta-md",
            help="jobname on the pbs server")
    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")
    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()

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

    params["MD.FinalTimeStep"] = args.mdstep
    params["MD.LengthTimeStep"] = args.timestep
    params["MD.InitialTemperature"] = args.initial_temp
    params["MD.TargetTemperature"] = args.target_temp
    params["MD.VariableCell"] = args.vc
    
    task = md_run()
    task.get_xyz(args.file)
    task.set_params(params=params)
    task.set_kpoints(kpoints_mp=args.kpoints_mp)
    task.md(directory=args.directory, runopt=args.runopt, mpi=args.mpi)


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
            ctl.submit(workdir=args.directory, jobfile="molecular-dynamics.pbs", server="pbs")
        elif args.server == "yh":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_yh.conf"))
            ctl.login()
            ctl.submit(workdir=args.directory, jobfile="molecular-dynamics.sub", server="yh")
