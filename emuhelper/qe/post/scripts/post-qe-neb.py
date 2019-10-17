#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.post.neb import neb_post



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously neb running directory", type=str, default="tmp-qe-neb")
    parser.add_argument("--nebint", help="xxx.int", type=str, default="pwscf.int")
    parser.add_argument("--nebdat", help="xxx.dat", type=str, default="pwscf.dat")

    args = parser.parse_args()

    task = neb_post()
    task.min_energy_path_gp(directory=args.directory, nebint=args.nebint, nebdat=args.nebdat, runopt="genrun")
