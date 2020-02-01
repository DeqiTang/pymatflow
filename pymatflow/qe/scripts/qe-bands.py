#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.static import static_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync

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

    parser.add_argument("--kpoints-option", type=str, default="crystal_b", 
            choices=["automatic", "gamma", "crystal_b"],
            help="Kpoints generation scheme option for band calculation")

    parser.add_argument("--kpoints-mp", type=str, default="1 1 1 0 0 0",
            help="Monkhorst-Pack kpoint grid, in format like '1 1 1 0 0 0'")

    parser.add_argument("--crystal-b-from", type=str, default="seekpath",
            choices=["seekpath", "manual"],
            help="the crystal_b kpoints from ? canbe seekpath or manual")
    
    parser.add_argument("--crystal-b-manual", type=str, nargs="+", default=None,
            help="manual input crystal_b, like --crystal-b-manual '0.000000 0.000000 0.000000 5  GAMMA' '0.500000 0.000000 0.000000 5  X'")

    parser.add_argument("--crystal-b-manual-file", type=str, default='crystal-b.txt',
            help="manual input crystal_b read from the file")

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
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")
    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")
    parser.add_argument("--jobname", type=str, default="qe-bands-structure",
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
    # process tpiba_b_manual
    crystal_b_manual = []
    if args.crystal_b_from == 'manual' and args.crystal_b_manual != None:
       # crystal_b_manual read from script argument args.crystal_b_manual
       for point in args.crystal_b_manual:
           crystal_b_manual.append([
               float(point.split()[0]), 
               float(point.split()[1]), 
               float(point.split()[2]),
               int(point.split()[3]),
               point.split()[4].upper()
               ])
    elif args.crystal_b_from == 'manual' and args.crystal_b_manual == None:
        # crystal_b_manual read from a file contains manual K_POINTS crystal_b setting
        # crystal_b_manual_file in format like this:
        # K_POINTS crystal_b
        # 2
        # 0.000000 0.000000 0.000000 5  #GAMMA
        # 0.500000 0.000000 0.000000 5  #X
        # 0.000000 0.500000 0.000000 5  #Y
        with open(args.crystal_b_manual_file, 'r') as fin:
          fin.readline()
          fin.readline()
          for point in fin:
              crystal_b_manual.append([
                  float(point.split()[0]), 
                  float(point.split()[1]), 
                  float(point.split()[2]),
                  int(point.split()[3]),
                  point.split()[4].split("#")[1].upper()
               ])

    # --------------------------------------------------------------------
    task = static_run()
    task.get_xyz(args.file)
    task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp, crystal_b_from=args.crystal_b_from, crystal_b_manual=crystal_b_manual)
    task.set_params(control=control, system=system, electrons=electrons)
    task.set_bands(bands_input=bands)
    task.bands(directory=args.directory, mpi=args.mpi, runopt=args.runopt, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)

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
            ctl.submit(workdir=args.directory, jobfile="band-structure.pbs", server="pbs")
        elif args.server == "yh":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_yh.conf"))
            ctl.login()
            ctl.submit(workdir=args.directory, jobfile="band-structure.sub", server="yh")
