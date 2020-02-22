#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.abinit.post.md import md_post



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously molecular dynamics running directory", type=str, default="tmp-abinit-md")
    parser.add_argument("--mdout", help="output file of molecular dynamics calculation", type=str, default="molecular-dynamics.out")

    args = parser.parse_args()

    os.chdir(args.directory)
    task = md_post(args.mdout)
    task.export()
    os.chdir("../")
