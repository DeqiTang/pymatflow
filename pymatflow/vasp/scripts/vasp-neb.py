#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse



from pymatflow.vasp.neb import neb_run

"""
usage:
指定的核数必须要是--nimage也就是INCAR中IMAGES变量的整数倍
否则MPI会报错, 这是VTST要求的.
"""

params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-vasp-neb-vtst",
        help="directory of the neb running")

    parser.add_argument("--runopt", type=str, default="gen",
            help="gen, run, or genrun")

    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    parser.add_argument("--nimage", type=int, default=5,
            help="number of image to interpolate. total image will be nimage+2.")

    parser.add_argument("--images", type=str, nargs="+",
            help="the image xyz file(--images first.xyz final.xyz)")

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

    parser.add_argument("--ediffg", type=float, default=-0.01,
            help="EDIFFG, default value: 10*EDIFF")

    parser.add_argument("--ibrion", type=int, default=3,
            choices=[1, 2, 3],
            help="IBRION = 1(ionic relaxation:RMM-DIIS[Quisi-Newton]), 2(ionic relaxation:CG), 3(damped molecular dynamics): refer to https://cms.mpi.univie.ac.at/wiki/index.php/IBRION for how to set the algorithm of optimization you need!")

    parser.add_argument("--isif", type=int, default=2,
            choices=[0, 1, 2, 3, 4, 5, 6, 7],
            help="ISIF = 0-7: refer to https://cms.mpi.univie.ac.at/wiki/index.php/ISIF for how to set the type of Geometri Optimization you need!")

    parser.add_argument("--potim", type=float, default=0,
            help="step width scaling (ionic relaxations), default: 0 -> use VTST optimizer")

    # --------------------------------------------------------------------------
    #                     neb related PARAMETERS
    # --------------------------------------------------------------------------
    parser.add_argument("--iopt", type=int, default=1,
            choices=[0, 1, 2],
            help="chioce for optimizer: 0->vasp, 1, 2->vtst")

    parser.add_argument("--lclimb", type=str, default="T",
            choices=["T", "F"],
            help="whether use climbing image")

    parser.add_argument("--spring", type=int, default=-5,
            help="gives the spring constant between the images as used in the elastic band method")


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

    parser.add_argument("--jobname", type=str, default="vasp-neb",
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

    params["EDIFFG"] = args.ediffg
    params["NSW"] = args.nsw
    params["IBRION"] = args.ibrion
    params["ISIF"] = args.isif
    params["POTIM"] = args.potim

    params["IOPT"] = args.iopt
    params["LCLIMB"] = args.lclimb
    params["SPRING"] = args.spring
    params["IMAGES"] = args.nimage

    task = neb_run()
    task.get_images(args.images)
    task.set_params(params=params)
    task.set_kpoints(kpoints_mp=args.kpoints_mp)
    task.nimage = args.nimage
    #task.neb(directory=args.directory, runopt=args.runopt, mpi=args.mpi, electrons=electrons_params, ions=ions_params, properties=properties_params, kpoints_mp=kpoints_mp)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.neb(directory=args.directory, runopt=args.runopt, auto=args.auto)
