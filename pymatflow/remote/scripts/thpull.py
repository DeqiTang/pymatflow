#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import sys

import argparse
from pymatflow.remote.ssh import Ssh
from pymatflow.remote.rsync import Rsync


"""
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", type=str,
            help="directory to pull from origin to local")

    parser.add_argument("--exclude", type=str, nargs="+",
            default=[],
            help="files or dirs to exclude from rsyncing. in format like this --excude '/xxx1' '/xxx2', (wildcard is supported!), default is copying everything")

    args = parser.parse_args()
    # server handle

    mover = Rsync()
    mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_yh.conf"))

    mover.copy(source=os.path.join(mover.user+"@"+mover.ip+":"+mover.serverdir, args.directory)+"/", target=args.directory, exclude=args.exclude)
