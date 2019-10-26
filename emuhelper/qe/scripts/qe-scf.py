#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.static import static_run

"""
usage:
    qe-scf.py -f xxx.xyz -k '2 2 2 0 0 0' --ecutwfc 100
"""


control_params = {}
system_params = {}
electrons_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of the static running", type=str, default="tmp-qe-static")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--runopt", help="gen, run, or genrun", type=str, default="genrun")
    parser.add_argument("--mpi", help="MPI command", type=str, default="")
    parser.add_argument("--ecutwfc", help="ecutwfc, default value: 100 Ry", type=int, default=100)
    parser.add_argument("--ecutrho", help="ecutrho, default value: 400 Ry", type=int, default=400)
    parser.add_argument("--kpoints-option", help="kpoints option", type=str, default="automatic")  
    parser.add_argument("-k", "--kpoints", help="set kpoints like '3 3 3 0 0 0'", type=str, default="3 3 3 0 0 0")
    parser.add_argument("--conv-thr", help="conv_thr", type=float, default=1.0e-6)
    parser.add_argument("--occupations", help="occupation type", type=str, default="smearing")
    parser.add_argument("--smearing", help="smearing type(gaussian, methfessel-paxton, mazari-vanderbilt, fermi-dirac), default: gaussian", type=str, default="gaussian")
    parser.add_argument("--degauss", help="value of the gaussian spreading (Ry) for brillouin-zone integration in metals.(defualt: 0.001 Ry)", type=float, default=0.001)
    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
    # ==========================================================   
    args = parser.parse_args()
    xyzfile = args.file
    system_params["ecutwfc"] = args.ecutwfc
    system_params["ecutrho"] = args.ecutrho
    system_params["occupations"] = args.occupations
    system_params["smearing"] = args.smearing
    system_params["degauss"] = args.degauss
    electrons_params["conv_thr"] = args.conv_thr
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]
    

    task = static_run(xyzfile)
    task.scf(directory=args.directory, runopt=args.runopt, mpi=args.mpi, control=control_params, system=system_params, electrons=electrons_params, kpoints_option=args.kpoints_option, kpoints_mp=kpoints_mp)
