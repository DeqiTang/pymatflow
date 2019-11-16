#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from pymatflow.qe.post.pdos import pdos_post

"""
usage: post-qe-pdos.py 
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously projwfc running directory", type=str, default="tmp-qe-static")
    parser.add_argument("--filpdos", help="filpdos name used in projwfc running", type=str, default="projwfc")

    args = parser.parse_args()

    directory = args.directory
    filpdos = args.filpdos

    task = pdos_post()
    task.get_data(directory=directory, filpdos=filpdos)
    task.export(directory=directory)
