#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
from pymatflow.qe.static import static_run
from pymatflow.remote.server import server_handle

"""
usage:
    qe-ir-raman.py xxx.xyz
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # for server
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.emuhelper/server.conf")

    args = parser.parse_args()

    task = static_run(sys.argv[1])
    task.ir_raman()

