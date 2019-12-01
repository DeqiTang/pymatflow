#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import argparse

from pymatflow.abinit.md import md_run

"""
usage: abinit-md.py xxx.xyz
"""

electrons_params = {}
kpoints_params = {}
ions_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-abinit-md",
            help="Directory for the molecular dynamics running.")

    parser.add_argument("-f", "--file", type=str,
            help="The xyz file name.")

    parser.add_argument("--runopt", type=str, default="genrun",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    # -----------------------------------------------------------------
    #                        scf related parameters
    # -----------------------------------------------------------------
    parser.add_argument("--ecut", type=int, default=15,
            help="Kinetic energy cutoff for wave functions in unit of Hartree, default value: 15 Hartree. for more information, refer to https://docs.abinit.org/variables/basic/#ecut")

    parser.add_argument("--ixc", type=int, default=11,
            choices=[1, 2, 3 ,4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 40, 41, 42],
            help="type of exchage-correlation functional. for more information, refer to https://docs.abinit.org/variables/basic/#ixc")
 
    parser.add_argument("--kptopt", type=int, default=1,
            help="Kpoints Generation scheme option: 0, 1, 2, 3, 4 or a negative value. for more information, refer to https://docs.abinit.org/variables/basic/#kptopt")

    parser.add_argument("--ngkpt", nargs="+", type=int,
            default=[1, 1, 1],
            help="number of grid points for kpoints generation. for more information, refer to https://docs.abinit.org/variables/basic/#ngkpt")

    parser.add_argument("--vdw-xc", type=int,
            default=None,
            choices=[0, 1, 2, 5, 6, 7, 10, 11, 14],
            help="Van Der Waals exchange-correlation functional. 5: DFT-D2, 6: DFT-D3, 7: DFT-D3(BJ). for more information, refer to https://docs.abinit.org/variables/vdw/#vdw_xc")

    parser.add_argument("--vdw-tol", type=float,
            default=None,
            help="Van Der Waals tolerance, only work when vdw_xc == 5 or 6 or 7. to be included in the potential a pair of atom must have contribution to the energy larger than vdw_tol. default value is 1.0e-10. for more information, refer to https://docs.abinit.org/variables/vdw/#vdw_tol")

    # -----------------------------------------------------------
    #                        ions moving related parameters
    # -----------------------------------------------------------
    parser.add_argument("--ionmov", type=int, default=6,
            choices=[1, 6, 7, 8, 9, 12, 13, 14, 24, 25],
            help="type of molecular dynamics algorithm, can be: 1 6 7 8 9 12 13 14 24 25. for more information, refer to https://docs.abinit.org/variables/rlx/#ionmov")

    parser.add_argument("--dtion", type=float, default=100,
            help="delta time for ions, in atom unit(One atomic time unit is 2.418884e-17 seconds), default=100. for more information, refer to https://docs.abinit.org/variables/rlx/#dtion.")

    parser.add_argument("--ntime", type=int, default=1000,
            help="number of time steps. for more information, refer to https://docs.abinit.org/variables/rlx/#ntime.")

    parser.add_argument("--mdtemp", type=float, nargs="+", default=[300, 300],
            help="molecular dynamics temperature. for more information, refer to https://docs.abinit.org/variables/rlx/#mdtemp.")

    parser.add_argument("--optcell", type=int,
            default=0,
            help="whether to optimize the cell shape and dimension")

    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
    # ==========================================================   
    args = parser.parse_args()

    electrons_params["ecut"] = args.ecut
    electrons_params["ixc"] = args.ixc
    electrons_params["vdw_xc"] = args.vdw_xc
    electrons_params["vdw_tol"] = args.vdw_tol

    kpoints_params["kptopt"] = args.kptopt
    kpoints_params["ngkpt"] = args.ngkpt

    ions_params["ionmov"] = args.ionmov
    ions_params["dtion"] = args.dtion
    ions_params["ntime"] = args.ntime
    ions_params["mdtemp(1)"] = args.mdtemp[0]
    ions_params["mdtemp(2)"] = args.mdtemp[1]
    ions_params["optcell"] = args.optcell

    task = md_run(args.file)
    task.md(directory=args.directory, mpi=args.mpi, runopt=args.runopt, electrons=electrons_params, kpoints=kpoints_params, ions=ions_params)

