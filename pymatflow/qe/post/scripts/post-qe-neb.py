#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.post.neb import neb_post



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-neb",
            help="previously neb running directory")

    parser.add_argument("--nebint", type=str, default="pwscf.int",
            help="xxx.int")

    parser.add_argument("--nebdat", type=str, default="pwscf.dat",
            help="xxx.dat")

    parser.add_argument("--inpname", type=str, default="min-energy-path.gp",
            help="inpname for the gnuplot script")

    parser.add_argument("--md", type=str, default="neb-report.md",
            help="Markdown report file name")

    parser.add_argument("--nebout", type=str, default="neb.out",
            help="output file of neb calculation")

    args = parser.parse_args()

    os.chdir(args.directory)
    task = neb_post(nebout=args.nebout)
    os.chdir("../")
    task.export(directory=args.directory, nebint=args.nebint, nebdat=args.nebdat, md=args.md)
