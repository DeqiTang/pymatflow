#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse


from pymatflow.vasp.md import md_run

"""
usage:
"""


params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", help="directory of the md running", type=str, default="tmp-vasp-md")

    parser.add_argument("-f", "--file", help="the xyz file name", type=str)

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

    parser.add_argument("--kpoints-mp", type=int, nargs="+",
            default=[1, 1, 1, 0, 0, 0],
            help="set kpoints like -k 1 1 1 0 0 0")

    parser.add_argument("--ismear", type=int, default=0,
            help="smearing type(methfessel-paxton(>0), gaussian(0), fermi-dirac(-1), tetra(-4), tetra-bloch-dorrected(-5)), default: 0")

    parser.add_argument("--sigma", type=float, default=0.01,
            help="determines the width of the smearing in eV.")

    parser.add_argument("--ivdw", type=int, default=None,
            choices=[0, 11, 12, 21, 202, 4],
            help="IVDW = 0(no correction), 1(dft-d2), 11(dft-d3 Grimme), 12(dft-d3 Becke-Jonson), default: None which means 0, no correction")
    # -----------------------------------------------------------------

    # -------------------------------
    # Ensemble and Thermostat
    # -------------------------------
    parser.add_argument("--ensemble", type=int, default=0,
            choices=[0, 1, 2, 3],
            help="choose ensemble: 0-NVE, 1-NVT, 2-NPT, 3-NPH")

    parser.add_argument("--thermostat", type=int, default=0,
            choices=[0, 1, 2, 3],
            help="choose thermostat: 0-Andersen, 1-Nose-Hoover, 2-Langevin, 3-Multiple Andersen.")
    # ----------------------
    # properties parametes
    # ---------------------
    #parser.add_argument("--lorbit", help="", type=int, default=None)
    #parser.add_argument("--loptics", help="", type=str, default="FALSE")


    # -----------------------------------------------------------------
    #                      run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh", "lsf_sz"],
            help="type of remote server, can be pbs or yh or lsf_sz")

    parser.add_argument("--jobname", type=str, default="vasp-md",
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
    params["ENCUT"] = args.encut
    params["EDIFF"] = args.ediff


    task = md_run()
    task.get_xyz(args.file)
    task.set_params(params=params)
    task.set_kpoints(kpoints_mp=args.kpoints_mp)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.md(directory=args.directory, runopt=args.runopt, auto=args.auto, ensemble=args.ensemble, thermostat=args.thermostat)
