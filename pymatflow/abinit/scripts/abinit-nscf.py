#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import argparse
from pymatflow.abinit.static import static_run

"""
usage:
"""

electrons_params = {}
kpoints_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-abinit-static",
            help="Directory for the static running.")

    parser.add_argument("-f", "--file", type=str,
            help="The xyz file name, containing the structure to be simulated")

    parser.add_argument("--runopt", type=str, default="genrun",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--properties", nargs="+", type=int,
            default=[],
            help="options for properties calculation")

    parser.add_argument("--iscf", type=int, default=-3,
            choices=[-3, -2, -1],
            help="set scf or nscf type. for more information, refer to https://docs.abinit.org/variables/basic/#iscf")

    parser.add_argument("--ecut", type=int, default=15,
            help="Kinetic energy cutoff for wave functions in unit of Hartree, default value: 15 Hartree. for more information, refer to https://docs.abinit.org/variables/basic/#ecut")

    parser.add_argument("--ixc", type=int, default=11,
            choices=[1, 2, 3 ,4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 40, 41, 42],
            help="type of exchage-correlation functional. for more information, refer to https://docs.abinit.org/variables/basic/#ixc")

    parser.add_argument("--kptopt", type=int, default=1,
            help="Kpoints Generation scheme option: 0, 1, 2, 3, 4 or a negative value. for more information, refer to https://docs.abinit.org/variables/basic/#kptopt")

    parser.add_argument("--ngkpt", nargs="+", type=int, default=[3, 3, 3],
            help="number of grid points for kpoints generation. for more information, refer to https://docs.abinit.org/variables/basic/#ngkpt")

    
    parser.add_argument("--prtdos", type=int, default=None,
            choices=[1, 2],
            help="prtdos. for more information, refer to https://docs.abinit.org/variables/files/#prtdos")

    parser.add_argument("--vdw-xc", type=int,
            default=None,
            choices=[0, 1, 2, 5, 6, 7, 10, 11, 14],
            help="Van Der Waals exchange-correlation functional. 5: DFT-D2, 6: DFT-D3, 7: DFT-D3(BJ). for more information, refer to https://docs.abinit.org/variables/vdw/#vdw_xc")

    parser.add_argument("--vdw-tol", type=float,
            default=None,
            help="Van Der Waals tolerance, only work when vdw_xc == 5 or 6 or 7. to be included in the potential a pair of atom must have contribution to the energy larger than vdw_tol. default value is 1.0e-10. fore more information, refer to https://docs.abinit.org/variables/vdw/#vdw_tol")

    
    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
    # ==========================================================   
    args = parser.parse_args()

    electrons_params["ecut"] = args.ecut
    electrons_params["ixc"] = args.ixc
    electrons_params["vdw_xc"] = args.vdw_xc
    electrons_params["vdw_tol"] = args.vdw_tol
    electrons_params["iscf"] = args.iscf

    kpoints_params["kptopt"] = args.kptopt
    kpoints_params["ngkpt"] = args.ngkpt

    electrons_params["prtdos"] = args.prtdos
    
    task = static_run(args.file)
    task.nscf(directory=args.directory, mpi=args.mpi, runopt=args.runopt, electrons=electrons_params, kpoints=kpoints_params, properties=args.properties)