#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.abinit.post.neb import neb_post



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously neb running directory", type=str, default="tmp-abinit-neb")
    parser.add_argument("--nebout", help="output file of neb calculation", type=str, default="neb.out")

    args = parser.parse_args()

    os.chdir(args.directory)
    task = neb_post(args.nebout)
    task.export()
    os.chdir("../")
