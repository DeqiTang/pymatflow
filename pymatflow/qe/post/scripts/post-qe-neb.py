#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.post.neb import neb_post



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously neb running directory", type=str, default="tmp-qe-neb")
    parser.add_argument("--nebint", help="xxx.int", type=str, default="pwscf.int")
    parser.add_argument("--nebdat", help="xxx.dat", type=str, default="pwscf.dat")
    parser.add_argument("--inpname", help="inpname for the gnuplot script", type=str, default="min-energy-path.gp")
    parser.add_argument("--md", help="Markdown report file name", type=str, default="neb-report.md")
    parser.add_argument("--nebout", help="output file of neb calculation", type=str, default="neb.out")

    args = parser.parse_args()

    task = neb_post()
    task.export(directory=args.directory, nebint=args.nebint, nebdat=args.nebdat, md=args.md, nebout=args.nebout)
