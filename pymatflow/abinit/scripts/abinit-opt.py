#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import argparse

from pymatflow.abinit.opt import opt_run

"""
usage:
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-abinit-opt",
            help="Directory to do the optimization calculation")

    parser.add_argument("-f", "--file", type=str,
            help="The xyz structure file with second line specifying cell parameters")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    # --------------------------------------------------------------------------
    parser.add_argument("--chkprim", type=int, default=1,
            choices=[0, 1],
            help="check whether the input cell is primitive. if your cell is not primitive, set chkprim to 0. for more information, refer to https://docs.abinit.org/variables/gstate/#chkprim")



    parser.add_argument("--ecut", type=int, default=15,
            help="Kinetic energy cutoff for wave functions in unit of Hartree, default value: 15 Hartree. for more information refer to https://docs.abinit.org/variables/basic/#ecut")

    parser.add_argument("--ixc", type=int, default=11,
            choices=[1, 2, 3 ,4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 40, 41, 42],
            help="type of exchage-correlation functional. fore more information, refer to https://docs.abinit.org/variables/basic/#ixc for more information")

    parser.add_argument("--kptopt", type=int, default=1,
            help="Kpoints Generation scheme option: 0, 1, 2, 3, 4 or a negative value. for more information, refer to https://docs.abinit.org/variables/basic/#kptopt")

    parser.add_argument("--ngkpt", nargs="+", type=int,
            default=[1, 1, 1],
            help="number of grid points for kpoints generation. for more information, refer to https://docs.abinit.org/variables/basic/#ngkpt")

    # vdw related parameters
    parser.add_argument("--vdw-xc", type=int, default=None,
            choices=[0, 1, 2, 5, 6, 7, 10, 11, 14],
            help="Van Der Waals exchange-correlation functional. 0: no correction, 1: vdW-DF1, 2: vdW-DF2, 5: DFT-D2, 6: DFT-D3, 7: DFT-D3(BJ). for more information, refer to https://docs.abinit.org/variables/vdw/#vdw_xc")

    parser.add_argument("--vdw-tol", type=float,
            default=None,
            help="Van Der Waals tolerance, only work when vdw_xc == 5 or 6 or 7. to be included in the potential a pair of atom must have contribution to the energy larger than vdw_tol. default value is 1.0e-10. fore more information, refer to https://docs.abinit.org/variables/vdw/#vdw_tol")

    # -----------------------------------------------------------
    #                        ions moving related parameters
    # -----------------------------------------------------------

    parser.add_argument("--ionmov", type=int, default=3,
            choices=[2, 3, 4, 5],
            help="type of ionmov algorithm. fore more information, refer to https://docs.abinit.org/variables/rlx/#ionmov")

    parser.add_argument("--optcell", type=int,
            choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            default=0,
            help="whether to optimize the cell shape and dimension. fore more information, refer to https://docs.abinit.org/variables/rlx/#optcell")

    parser.add_argument("--ecutsm", type=float, default=None,
            help="when optcell != 0, must specify encutsm larser than zero. for more information refer to https://docs.abinit.org/variables/rlx/#ecutsm")


    # -----------------------------------------------------------------
    #                      run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="opt-cubic",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
    # ==========================================================
    args = parser.parse_args()

    params = {}
    kpoints = {}

    params["chkprim"] = args.chkprim
    params["ecut"] = args.ecut
    params["ixc"] = args.ixc
    params["vdw_xc"] = args.vdw_xc
    params["vdw_tol"] = args.vdw_tol

    kpoints["kptopt"] = args.kptopt
    kpoints["ngkpt"] = args.ngkpt

    params["ionmov"] = args.ionmov
    params["optcell"] = args.optcell
    params["ecutsm"] = args.ecutsm

    task = opt_run()
    task.get_xyz(args.file)
    task.set_params(params=params)
    task.set_kpoints(kpoints=kpoints)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.optimize(directory=args.directory, runopt=args.runopt, auto=args.auto)
