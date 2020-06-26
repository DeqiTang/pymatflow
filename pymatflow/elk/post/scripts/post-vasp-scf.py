#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.vasp.post.scf import scf_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of relax running", type=str, default="tmp-vasp-static")

    args = parser.parse_args()

    os.chdir(args.directory)
    task = scf_post()
    task.get_outcar("OUTCAR")
    task.export()
