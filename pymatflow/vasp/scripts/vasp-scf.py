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

    parser.add_argument("--lorbit", type=int, default=None,
            choices=[0, 1, 2, 5, 10, 11, 12],
            help="together with an appropriate RWIGS, determines whether the PROCAR or PROOUT files are written")

    parser.add_argument("--loptics", type=str, default="FALSE",
            choices=["TRUE", "FALSE"],
            help="calculates the frequency dependent dielectric matrix after the electronic ground state has been determined.")

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
    # special
    parser.add_argument("--algo", type=str, default=None,
            choices=["N", "D", "V", "F"],  #"Exact", "G0W0", "GW0", "GW"],
            help=" a convenient option to specify the electronic minimisation algorithm (as of VASP.4.5) and/or to select the type of GW calculations")

    # hybrid functional
    parser.add_argument("--lhfcalc", type=str, default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help=" specifies whether Hartree-Fock/DFT hybrid functional type calculations are performed")

    parser.add_argument("--hfscreen", type=float, default=None,
            choices=[0.3, 0.2],
            help=" specifies the range-separation parameter in range separated hybrid functionals: HSE03->0.3, HSE06->0.2, must also set LHFCALC=.TRUE.")

    parser.add_argument("--lsubrot", type=str, default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help="This flag can be set for hybrid functionals (HF-type calculations).")

    # properties parameters
    parser.add_argument("--lelf", type=str, default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help="LELF determines whether to create an ELFCAR file or not.")

    # -----------------------------------------------------------------
    #                      run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh", "lsf_sz"],
            help="type of remote server, can be pbs or yh or lsf_sz")

    parser.add_argument("--jobname", type=str, default="vasp-scf",
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
    params["LORBIT"] = args.lorbit
    params["LOPTICS"] = args.loptics
    params["LSUBROT"] = args.lsubrot


    params["ISPIN"] = args.ispin
    params["MAGMOM"] = args.magmom # magmom can be a list that can be automatically dealt with by base.incar.to_incar()
    params["LNONCOLLLINEAR"] = args.lnoncollinear
    params["LSORBIT"] = args.lsorbit
    params["ALGO"] = args.algo
    params["LHFCALC"] = args.lhfcalc
    params["HFSCREEN"] = args.hfscreen

    params["LELF"] = args.lelf

    task = static_run()
    task.get_xyz(args.file)
    task.set_params(params)
    task.set_kpoints(kpoints_mp=args.kpoints_mp)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto)
