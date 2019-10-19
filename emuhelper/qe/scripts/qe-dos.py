#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.static import static_run

"""
usage:
    qe-single-point.py -f xxx.xyz
Note:
    do not use --bzsum option now, because there
    is some problem with static_run.dos() when 
    dealring with bz_sum in input file
    namelist bug in dos.x when bz_sum is present
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)
    parser.add_argument("-d", "--directory", help="directory for the calculation", type=str, default="tmp-qe-static")
    parser.add_argument("--fildos", help="output dos file name", type=str, default="dosx.dos")
    parser.add_argument("--bzsum", help="brillouin summation type", type=str, default="smearing")
    parser.add_argument("--ngauss", help="gaussian broadening type", type=str, default='default')
    parser.add_argument("--degauss", help="gaussian broadening", type=str, default='default')
    parser.add_argument("--emin", help="min energy for DOS", type=str, default='default')
    parser.add_argument("--emax", help="max energy for DOS", type=str, default='default')
    parser.add_argument("--deltae", help="DeltaE: energy grid step (eV)", type=str, default='default')

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file
    directory = args.directory
    fildos = args.fildos
    bzsum = args.bzsum
    if args.ngauss == 'default':
        ngauss = args.ngauss
    else:
        ngauss = int(args.ngauss)
    if args.degauss == 'default':
        degauss = args.degauss
    else:
        degauss = float(args.degauss)
    if args.emin == 'default':
        emin = args.emin
    else:
        emin = float(args.emin)
    if args.emax == 'default':
        emax = args.emax
    else:
        emax = float(args.emax)
    if args.deltae == 'default':
        deltae = args.deltae
    else:
        deltae = float(args.deltae)


    task = static_run(xyzfile)
    task.dos(directory=directory, fildos=fildos, bz_sum=bzsum, ngauss=ngauss, degauss=degauss, emin=emin, emax=emax, deltae=deltae)
