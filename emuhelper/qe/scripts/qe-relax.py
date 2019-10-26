#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.opt import opt_run

"""
usage: qe-relax.py -f xxx.xyz
"""

control_params = {}
system_params = {}
electrons_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory for the relax running", type=str, default="tmp-qe-relax")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--runopt", help="gen, run, or genrun", type=str, default="genrun")
    parser.add_argument("--mpi", help="MPI command", type=str, default="")
    parser.add_argument("--etot-conv-thr", type=float, default=1.0e-4)
    parser.add_argument("--forc-conv-thr", type=float, default=1.0e-3)
    parser.add_argument("--nstep", help="number of maximum optimization steps", type=int, default=50)
    parser.add_argument("--ecutwfc", help="ecutwfc", type=int, default=100)
    parser.add_argument("--ecutrho", help="ecutrho", type=int, default=400)
    parser.add_argument("--kpoints-option", help="kpoints option", type=str, default="automatic")
    parser.add_argument("-k", "--kpoints", help="set kpoints like '1 1 1 0 0 0'", type=str, default="1 1 1 0 0 0")
    parser.add_argument("--conv-thr", help="the conv_thr for scf", type=float, default=1.0e-6)
    parser.add_argument("--occupations", help="occupation type", type=str, default="smearing")
    parser.add_argument("--smearing", help="smearing type", type=str, default="gaussian")
    parser.add_argument("--degauss", help="value of the gaussian spreading (Ry) for brillouin-zone integration in metals.", type=float, default=0.001)
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file
    control_params["etot_conv_thr"] = args.etot_conv_thr
    control_params["forc_conv_thr"] = args.forc_conv_thr
    control_params["nstep"] = args.nstep
    system_params["ecutwfc"] = args.ecutwfc
    system_params["ecutrho"] = args.ecutrho
    system_params["occupations"] = args.occupations
    system_params["smearing"] = args.smearing
    system_params["degauss"] = args.degauss
    electrons_params["conv_thr"] = args.conv_thr
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]
 
    task = opt_run(xyzfile)
    task.relax(directory=args.directory, runopt=args.runopt, mpi=args.mpi, control=control_params, system=system_params, electrons=electrons_params, kpoints_option=args.kpoints_option, kpoints_mp=kpoints_mp)
