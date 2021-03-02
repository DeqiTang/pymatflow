#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse


from pymatflow.vasp.static import static_run

"""
usage:
"""

params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-vasp-static",
            help="directory of the static running")

    parser.add_argument("-f", "--file", type=str,
            help="the xyz file name")

    parser.add_argument("--runopt", type=str, default="gen",
            help="gen, run, or genrun")

    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    # --------------------------------------------------------
    #                   INCAR PARAMETERS
    # --------------------------------------------------------
    parser.add_argument("--prec", type=str, default="Normal",
            choices=["Normal", "Accurate", "A", "N"],
            help="PREC, default value: Normal")

    parser.add_argument("--encut", type=int, default=300,
            help="ENCUT, default value: 300 eV")

    parser.add_argument("--ediff", type=float, default=1.0e-4,
            help="EDIFF, default value: 1.0e-4")

    parser.add_argument("--kpath", type=str, nargs="+", default=None,
            help="set kpoints for band structure calculation manually")

    parser.add_argument("--kpath-file", type=str, default="kpath-from-seekpath.txt",
            help="set kpoints for band structure calculation manually from file")

    parser.add_argument("--kpath-intersections", type=int, default=15,
            help="intersection of the line mode kpoint for band calculation")

    parser.add_argument("--ismear", type=int, default=0,
            help="smearing type(methfessel-paxton(>0), gaussian(0), fermi-dirac(-1), tetra(-4), tetra-bloch-dorrected(-5)), default: 0")

    parser.add_argument("--sigma", type=float, default=0.01,
            help="determines the width of the smearing in eV.")

    parser.add_argument("--ivdw", type=int, default=None,
            choices=[0, 11, 12, 21, 202, 4],
            help="IVDW = 0(no correction), 1(dft-d2), 11(dft-d3 Grimme), 12(dft-d3 Becke-Jonson), default: None which means 0, no correction")

    # magnetic related
    parser.add_argument("--ispin", type=int, default=None,
            choices=[1, 2],
            help="specifies spin polarization: 1->no spin polarized, 2->spin polarized(collinear). combine SIPIN with MAGMOM to study collinear magnetism.")

    parser.add_argument("--magmom", type=float, nargs="+", default=None,
            help="Specifies the initial magnetic moment for each atom, if and only if ICHARG=2, or if ICHARG=1 and the CHGCAR file contains no magnetisation density.")

    parser.add_argument("--lnoncollinear", type=str, default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help="specifies whether fully non-collinear magnetic calculations are performed")

    parser.add_argument("--lsorbit", type=str, default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help="specifies whether spin-orbit coupling is taken into account.")

    # -----------------------------------------------------------------


    # -----------------------------------------------------------------
    #                      run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh", "lsf_sz"],
            help="type of remote server, can be pbs or yh or lsf_sz")

    parser.add_argument("--jobname", type=str, default="vasp-bands",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
    # ==========================================================
    args = parser.parse_args()

    params["PREC"] = args.prec
    params["ENCUT"] = args.encut
    params["EDIFF"] = args.ediff
    params["ISMEAR"] = args.ismear
    params["SIGMA"] = args.sigma
    params["IVDW"] = args.ivdw

    params["ISPIN"] = args.ispin
    params["MAGMOM"] = args.magmom # magmom can be a list that can be automatically dealt with by base.incar.to_incar()
    params["LNONCOLLLINEAR"] = args.lnoncollinear
    params["LSORBIT"] = args.lsorbit

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


    task = static_run()
    task.get_xyz(args.file)
    task.set_params(params)
    # must set option to 'bands'
    task.set_kpoints(option="bands", kpath=kpath, kpath_intersections=args.kpath_intersections)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.bands(directory=args.directory, runopt=args.runopt, auto=args.auto)
