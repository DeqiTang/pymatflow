#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.static import static_run

"""
usage:
    qe-pdos.py -f xxx.xyz
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()    
    parser.add_argument("-d", "--directory", help="directory for the calculation", type=str, default="tmp-qe-static")
    parser.add_argument("-f", "--file", help="the xyz file", type=str)
    parser.add_argument("--filpdos", help="output projected dos file name", type=str, default="projwfc")
    parser.add_argument("--ngauss", help="gaussian broadening type", type=str, default='default')
    parser.add_argument("--degauss", help="gaussian broadening", type=str, default='default')
    parser.add_argument("--emin", help="min energy for DOS", type=str, default='default')
    parser.add_argument("--emax", help="max energy for DOS", type=str, default='default')
    parser.add_argument("--deltae", help="DeltaE: energy grid step (eV)", type=str, default='default')
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    directory = args.directory
    xyzfile = args.file
    filpdos = args.filpdos

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
    task.projwfc(directory=directory, filpdos=filpdos, ngauss=ngauss, degauss=degauss, emin=emin, emax=emax, deltae=deltae)
