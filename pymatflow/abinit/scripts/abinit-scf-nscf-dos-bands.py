#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import argparse

from pymatflow.abinit.static import static_run

"""
usage:
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-abinit-static",
            help="Directory to do the static scf calculation")

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

    parser.add_argument("--properties", nargs="+", type=int,
            default=[],
            help="options for properties calculation")

    parser.add_argument("--iscf", type=int, default=7,
            choices=[0, 1, 2, 3, 4, 5, 7, 12, 13, 14, 15, 17],
            help="set scf or nscf type. for more information, refer to https://docs.abinit.org/variables/basic/#iscf")

    parser.add_argument("--ecut", type=int, default=15,
            help="Kinetic energy cutoff for wave functions in unit of Hartree, default value: 15 Hartree. for more information, refer to https://docs.abinit.org/variables/basic/#ecut")

    parser.add_argument("--ixc", type=int, default=11,
            choices=[1, 2, 3 ,4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 40, 41, 42],
            help="type of exchage-correlation functional. for more information, refer to https://docs.abinit.org/variables/basic/#ixc")

    parser.add_argument("--kptopt", type=int, default=1,
            choices=[1],
            help="Kpoints Generation scheme option: 0, 1, 2, 3, 4 or a negative value. for more information, refer to https://docs.abinit.org/variables/basic/#kptopt")

    parser.add_argument("--ngkpt", nargs="+", type=int,
            default=[1, 1, 1],
            help="number of grid points for kpoints generation. for more information, refer to https://docs.abinit.org/variables/basic/#ngkpt")


    parser.add_argument("--kptbounds", type=str, nargs="+", default=None,
            help="manual input kpath for band structure calculation")

    parser.add_argument("--kptbounds-file", type=str, default="kpath-from-seekpath.txt",
            help="file to read the kpath for band structure calculation")

    # electron occupation
    parser.add_argument("--occopt", type=int, default=3,
            choices=[0, 1, 2, 3, 4, 5, 6, 7],
            help="Controls how input parameters nband, occ, and wtk are handled. for more information, refer to https://docs.abinit.org/variables/basic/#occopt")

    parser.add_argument("--nband", type=int, nargs="+", default=None,
            help="Gives number of bands, occupied plus possibly unoccupied, for which wavefunctions are being computed along with eigenvalues. for more information, refer to https://docs.abinit.org/variables/basic/#nband")

    parser.add_argument("--occ", nargs="+", type=float, default=None,
            help="Gives occupation numbers for all bands in the problem. Needed if occopt == 0 or occopt == 2. Ignored otherwise. Also ignored when iscf = -2. refer to https://docs.abinit.org/variables/gstate/#occ")

    # magnetic related parameters
    parser.add_argument("--nsppol", type=int, default=None,
            choices=[1, 2],
            help="Give the number of INDEPENDENT spin polarisations, for which there are non- related wavefunctions. Can take the values 1 or 2. for more information, refer to https://docs.abinit.org/variables/basic/#nsppol")

    # vdw related parameters
    parser.add_argument("--vdw-xc", type=int, default=None,
            choices=[0, 1, 2, 5, 6, 7, 10, 11, 14],
            help="Van Der Waals exchange-correlation functional. 0: no correction, 1: vdW-DF1, 2: vdW-DF2, 5: DFT-D2, 6: DFT-D3, 7: DFT-D3(BJ). for more information, refer to https://docs.abinit.org/variables/vdw/#vdw_xc")

    parser.add_argument("--vdw-tol", type=float,
            default=None,
            help="Van Der Waals tolerance, only work when vdw_xc == 5 or 6 or 7. to be included in the potential a pair of atom must have contribution to the energy larger than vdw_tol. default value is 1.0e-10. fore more information, refer to https://docs.abinit.org/variables/vdw/#vdw_tol")

    parser.add_argument("--prtden", type=int ,default=1,
            choices=[0, 1],
            help="print the density. for more information, refer to https://docs.abinit.org/variables/files/#prtden")

    parser.add_argument("--prtdos", type=int, default=None,
            choices=[0, 1, 2, 3],
            help="can be 0, 1, 2, 3. for more information, refer to https://docs.abinit.org/variables/files/#prtdos")

    # -----------------------------------------------------------------
    #                      run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="abinit-scf",
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

    params["occopt"] = args.occopt
    params["nband"] = args.nband
    params["occ"] = args.occ

    params["nsppol"] = args.nsppol

    params["prtden"] = args.prtden
    params["prtdos"] = args.prtdos


    # --------------------------------------------------------------------------
    # process kptbounds
    # --------------------------------------------------------------------------
    if args.kptbounds != None:
        # kptbounds from script argument args.kptbounds.
        kptbounds = []
        for kpoint in args.kptbounds:
            if kpoint.split()[4] != "|":
                kptbounds.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    int(kpoint.split()[4]),
                    ])
            elif kpoint.split()[4] == "|":
                kptbounds.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    "|",
                    ])
    elif args.kptbounds == None:
        # kptbounds read from file specified by args.kptbounds_file
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
        kptbounds = []
        with open(args.kptbounds_file, 'r') as fin:
            kptbounds_file = fin.readlines()
        nk = int(kptbounds_file[0])
        for i in range(nk):
            if kptbounds_file[i+1].split("\n")[0].split()[4] != "|":
                kptbounds.append([
                    float(kptbounds_file[i+1].split()[0]),
                    float(kptbounds_file[i+1].split()[1]),
                    float(kptbounds_file[i+1].split()[2]),
                    kptbounds_file[i+1].split()[3].split("#")[1].upper(),
                    int(kptbounds_file[i+1].split()[4]),
                    ])
            elif kptbounds_file[i+1].split("\n")[0].split()[4] == "|":
                kptbounds.append([
                    float(kptbounds_file[i+1].split()[0]),
                    float(kptbounds_file[i+1].split()[1]),
                    float(kptbounds_file[i+1].split()[2]),
                    kptbounds_file[i+1].split()[3].split("#")[1].upper(),
                    '|',
                    ])
    else:
        pass
    # --------------------------------------------------------------------





    task = static_run()
    task.get_xyz(args.file)
    task.set_params(params=params)
    task.set_kpoints(kpoints=kpoints)
    task.set_properties(properties=args.properties)
    task.dataset[3].electrons.kpoints.set_band(kptbounds=kptbounds)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto)
