#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.md import md_run

"""
usage: qe-md.py xxx.xyz
"""

control_params = {}
system_params = {}
electrons_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory for the vc-md running", type=str, default="tmp-qe-vc-md")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--runopt", help="gen, run, or genrun", type=str, default="genrun")
    parser.add_argument("--mpi", help="MPI command", type=str, default="")
    parser.add_argument("--nstep", help="maximum ion steps", type=int, default=50)
    parser.add_argument("--ecutwfc", help="ecutwfc", type=int, default=100)
    parser.add_argument("--ecutrho", help="ecutrho", type=int, defualt=400)
    parser.add_argument("-k", "--kpoints", help="set kpoints like '1 1 1 0 0 0'", type=str, default="1 1 1 0 0 0")
    parser.add_argument("--conv-thr", help="conv_thr of scf", type=float, default=1.e-6)
    parser.add_argument("--occupations", help="occupation type", type=str, default="smearing")
    parser.add_argument("--smearing", help="smearing type", type=str, default="gaussian")
    parser.add_argument("--degauss", help="value of the gaussian spreading (Ry) for brillouin-zone integration in metals.", type=float, default=0.001)
    parser.add_argument("--vdw-corr", help="vdw_corr = dft-d, dft-d3, ts-vdw, xdm", type=str, default='none')

    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
    # ==========================================================   

    args = parser.parse_args()
    xyzfile = args.file
    control_params["nstep"] = args.nstep
    system_params["ecutwfc"] = args.ecutwfc
    system_params["ecutrho"] = args.ecutrho
    system_params["occupations"] = args.occupations
    system_params["smearing"] = args.smearing
    system_params["degauss"] = args.degauss
    system_params["vdw_corr"] = args.vdw_corr
    electrons_params["conv_thr"] = args.conv_thr
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]
 
    task = md_run(xyzfile)
    task.vc_md(directory=args.directory, runopt=args.runopt, mpi=args.mpi, control=control_params, system=system_params, electrons=electrons_params, kpoints_mp=kpoints_mp)
