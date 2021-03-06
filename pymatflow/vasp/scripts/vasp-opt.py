#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse



from pymatflow.vasp.opt import opt_run

"""
usage:
"""


params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of the static running", type=str, default="tmp-vasp-optimization")

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

    parser.add_argument("--nsw", type=int, default=50,
            help="NSW sets the maximum number of ionic steps")

    parser.add_argument("--ediffg", type=float, default=None,
            help="EDIFFG, default value: 10*EDIFF")

    parser.add_argument("--ibrion", type=int, default=2,
            choices=[1, 2, 3],
            help="IBRION = 1(ionic relaxation:RMM-DIIS[Quisi-Newton]), 2(ionic relaxation:CG), 3(damped molecular dynamics): refer to https://cms.mpi.univie.ac.at/wiki/index.php/IBRION for how to set the algorithm of optimization you need!")

    parser.add_argument("--isif", type=int, default=2,
            choices=[0, 1, 2, 3, 4, 5, 6, 7],
            help="ISIF = 0-7: refer to https://cms.mpi.univie.ac.at/wiki/index.php/ISIF for how to set the type of Geometri Optimization you need!")

    parser.add_argument("--potim", type=float, default=None,
            help="step width scaling (ionic relaxations), default: None = 0.5 in opt")


    # special
    parser.add_argument("--algo", type=str, default=None,
            choices=["N", "D", "V", "F"],  #"Exact", "G0W0", "GW0", "GW"],
            help=" a convenient option to specify the electronic minimisation algorithm (as of VASP.4.5) and/or to select the type of GW calculations")

    parser.add_argument("--ialgo", type=int, default=None,
            choices=[5, 6, 7, 8, 38, 44, 46, 48],
            help="IALGO selects the algorithm used to optimize the orbitals.Mind: We strongly urge the users to select the algorithms via ALGO. Algorithms other than those available via ALGO are subject to instabilities.")

    parser.add_argument("--addgrid", type=str, default=None,
            choices=[".TRUE.", ".FALSE.", "T", "F"],
            help="ADDGRID determines whether an additional support grid is used for the evaluation of the augmentation charges.")

    parser.add_argument("--isym", type=int, default=None,
            choices=[-1, 0, 1, 2, 3],
            help=" ISYM determines the way VASP treats symmetry.")

    parser.add_argument('--lreal', type=str, default=None,
            choices=["T", "F", ".TRUE.", ".FALSE.", "O", "On", "A", "Auto"],
            help="LREAL determines whether the projection operators are evaluated in real-space or in reciprocal space.")


    # -----------------------------------------------------------------
    #                      run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh", "lsf_sz"],
            help="type of remote server, can be pbs or yh or lsf_sz")

    parser.add_argument("--jobname", type=str, default="vasp-opt",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")

    # pymatflow.vasp inside PARAMETERS
    parser.add_argument("--selective-dynamics", type=str, default="False",
            choices=["True", "False", "T", "F"],
            help="whether use selective dyanmics")

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
    params["EDIFFG"] = args.ediffg
    params["NSW"] = args.nsw
    params["IBRION"] = args.ibrion
    params["ISIF"] = args.isif
    params["POTIM"] = args.potim

    params["ALGO"] = args.algo
    params["IALGO"] = args.ialgo
    params["ADDGRID"] = args.addgrid
    params["ISYM"] = args.isym
    params["LREAL"] = args.lreal

    task = opt_run()
    task.get_xyz(args.file)
    task.set_params(params=params)
    task.set_kpoints(kpoints_mp=args.kpoints_mp)
    task.poscar.selective_dynamics = True if args.selective_dynamics.upper()[0] == "T" else False
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.optimize(directory=args.directory, runopt=args.runopt, auto=args.auto)
