#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse


from emuhelper.cp2k.opt import opt_run

"""
Usage:
    cp2k-geo-opt.py -f xxx.xyz
    xxx.xyz is the input structure file

"""

force_eval = {}
motion = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of the calculation", type=str, default="tmp-cp2k-geo-opt")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--cutoff", help="CUTOFF, default value: 100 Ry", type=int, default=100)
    parser.add_argument("--rel-cutoff", help="REL_CUTOFF, default value: 60 Ry", type=int , default=60)
    parser.add_argument("--xc-functional", help="XC_FUNCTIONAL type", type=str, default="PBE")
    parser.add_argument("-k", "--kpoints", help="set kpoints like '3 3 3 0 0 0'", type=str, default="3 3 3 0 0 0")
    parser.add_argument("--added-mos", help="ADDED_MOS for SCF", type=int, default=0)
    parser.add_argument("--smear", type=bool, default=False,
            choices=[True, False],
            help="switch on off smearing for occupation")
    parser.add_argument("--smear-method", help="smearing method", type=str, default="FERMI_DIRAC")
    parser.add_argument("--electronic-temp", help="ELECTRON_TEMPERATURE for FERMI_DIRAC SMEAR", type=float, default=300)
    parser.add_argument("--window-size", help="Size of the energy window centred at the Fermi level for ENERGY_WINDOW type smearing", type=float, default=0)
    # motion/geo_opt related
    parser.add_argument("--optimizer", type=str, default="BFGS",
            help="optimization algorithm for geometry optimization: BFGS, CG, LBFGS")
    parser.add_argument("--max-iter", type=int, default=200,
            help="maximum number of geometry optimization steps.")
    parser.add_argument("--type", type=str, default="MINIMIZATION",
            help="specify which kind of geometry optimization to perform: MINIMIZATION(default), TRANSITION_STATE")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()
    directory = args.directory
    xyzfile = args.file
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]
    force_eval["DFT-MGRID-CUTOFF"] = args.cutoff
    force_eval["DFT-XC-XC_FUNCTIONAL"] = args.xc_functional
    force_eval["DFT-SCF-ADDED_MOS"] = args.added_mos
    force_eval["DFT-SCF-SMEAR"] = args.smear
    force_eval["DFT-SCF-SMEAR-SMEAR-METHOD"] = args.smear
    force_eval["DFT-SCF-SMEAR-ELECTRONIC_TEMPERATURE"] = args.electronic_temp
    force_eval["DFT-SCF-SMEAR-WINDOW_SIZE"] = args.window_size
    
    motion["GEO_OPT-MAX_ITER"] = args.max_iter
    motion["GEO_OPT-OPTIMIZAER"] = args.optimizer
    motion["GEO_OPT-TYPE"] = args.type

    task = opt_run(xyzfile)
    task.geo_opt(directory="tmp-cp2k-geo-opt", runopt="genrun", force_eval=force_eval, motion=motion)

