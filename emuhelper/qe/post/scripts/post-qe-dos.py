#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.post.dos import dos_post

"""
usage: post-qe-dos.py xxx.dos
"""

if __name__ == "__main__":
   
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously projwfc running directory", type=str, default="tmp-qe-static")
    parser.add_argument("--fildos", help="fildos name used in dos.x running", type=str, default="dosx.dos")

    args = parser.parse_args()

    directory = args.directory
    fildos = args.fildos

    task = dos_post()
    task.get_data(directory=directory, fildos=fildos)
    task.export(directory=directory)
