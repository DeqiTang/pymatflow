#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.static import static_run

"""
usage:
    qe-nscf.py -f xxx.xyz -k '4 4 4 0 0 0'
"""

system_params = {}
electrons_params = {}
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory for the calcualtion", type=str, default="tmp-qe-static")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--mpi", help="MPI command", type=str, default="")
    parser.add_argument("--ecutwfc", help="ecutwfc, default value: 100 Ry", type=int, default=100)
    parser.add_argument("-k", "--kpoints", help="set kpoints like '4 4 4 0 0 0'", type=str, default="4 4 4 0 0 0")
    parser.add_argument("--occupations", help="occupation type", type=str, default="smearing")
    parser.add_argument("--smearing", help="smearing type", type=str, default="gaussian")
    parser.add_argument("--degauss", help="value of the gaussian spreading (Ry) for brillouin-zone integration in metals.", type=float, default=0.001)
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()
    xyzfile = args.file
    directory = args.directory
    system_params["ecutwfc"] = args.ecutwfc
    system_params["occupations"] = args.occupations
    system_params["smearing"] = args.smearing
    system_params["degauss"] = args.degauss
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]

    task = static_run(xyzfile)
    task.nscf(directory=directory, runopt='genrun', mpi=args.mpi, system=system_params, electrons=electrons_params, kpoints_mp=kpoints_mp)
