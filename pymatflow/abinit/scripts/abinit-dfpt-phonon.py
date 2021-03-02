#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import argparse

from pymatflow.remote.server import server_handle
from pymatflow.abinit.dfpt import dfpt_phonon

"""
usage:
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-abinit-dfpt-phonon",
            help="Directory for the dfpt elastic running.")

    parser.add_argument("-f", "--file", type=str,
            help="The xyz structure file with second line specifying cell parameters")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    # -------------------------------------------------------------------
    parser.add_argument("--chkprim", type=int, default=1,
            choices=[0, 1],
            help="check whether the input cell is primitive. if your cell is not primitive, set chkprim to 0. for more information, refer to https://docs.abinit.org/variables/gstate/#chkprim")

    parser.add_argument("--properties", nargs="+", type=int,
            default=[],
            help="options for properties calculation")

    parser.add_argument("--ecut", type=int, default=15,
            help="Kinetic energy cutoff for wave functions in unit of Hartree, default value: 15 Hartree. for more information, refer to https://docs.abinit.org/variables/basic/#ecut")

    parser.add_argument("--ixc", type=int, default=11,
            choices=[1, 2, 3 ,4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 40, 41, 42],
            help="type of exchage-correlation functional. for more information, refer to https://docs.abinit.org/variables/basic/#ixc")

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


    parser.add_argument("--kpath", type=str, nargs="+", default=None,
            help="set kpoints for band structure calculation manually")

    parser.add_argument("--kpath-file", type=str, default="kpath-from-seekpath.txt",
            help="set kpoints for band structure calculation manually from file")


    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

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

    kpoints["ngkpt"] = args.ngkpt


    # if band structure is in the printout option get the kpath
    if args.kpath != None:
        # kpath from script argument args.kpath
        kpath = []
        for kpoint in args.kpath:
            if kpoint.split()[4] != "|":
                kpath.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    int(kpoint.split()[4]),
                    ])
            elif kpoint.split()[4] == "|":
                kpath.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    "|",
                    ])
    elif args.kpath == None:
        # kpath read from file specified by args.kpath_file
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
        kpath = []
        with open(args.kpath_file, 'r') as fin:
            kpath_file = fin.readlines()
        nk = int(kpath_file[0])
        for i in range(nk):
            if kpath_file[i+1].split("\n")[0].split()[4] != "|":
                kpath.append([
                    float(kpath_file[i+1].split()[0]),
                    float(kpath_file[i+1].split()[1]),
                    float(kpath_file[i+1].split()[2]),
                    kpath_file[i+1].split()[3].split("#")[1].upper(),
                    int(kpath_file[i+1].split()[4]),
                    ])
            elif kpath_file[i+1].split("\n")[0].split()[4] == "|":
                kpath.append([
                    float(kpath_file[i+1].split()[0]),
                    float(kpath_file[i+1].split()[1]),
                    float(kpath_file[i+1].split()[2]),
                    kpath_file[i+1].split()[3].split("#")[1].upper(),
                    '|',
                    ])
    else:
        pass

    # -------------------------------------------------------------------

    task = dfpt_phonon()
    task.get_qpath(kpath)

    task.get_xyz(args.file)
    task.set_params(params=params)
    task.set_kpoints(kpoints=kpoints)
    task.set_properties(properties=args.properties)
    task.run(directory=args.directory, runopt=args.runopt, auto=args.auto)
