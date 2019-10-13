#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.static import static_run

"""
usage:
    qe-nscf.py -f xxx.xyz -k '4 4 4 0 0 0'
"""

system_params = {
        "ecutwfc": 150
        }

    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("--ecutwfc", help="ecutwfc", type=int, default=100)
    parser.add_argument("-k", "--kpoints", help="set kpoints like '4 4 4 0 0 0'", type=str, default="4 4 4 0 0 0")
    
    args = parser.parse_args()
    xyzfile = args.file
    system_params["ecutwfc"] = args.ecutwfc
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]

    task = static_run(xyzfile)
    task.nscf(runopt='genrun', system=system_params, kpoints_mp=kpoints_mp)
