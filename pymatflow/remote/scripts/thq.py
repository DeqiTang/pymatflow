#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os

import argparse
from pymatflow.remote.ssh import Ssh


"""
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", type=str, default="no",
            choices=["yes", "no"],
            help="check all the queue for all user on server. yes: check all, 0: only check yours")
    args = parser.parse_args()
    # server handle
    ctl = Ssh()
    ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_yh.conf"))
    ctl.login()
    if args.all == "yes":
        ctl.execute("yhqueue")
    elif args.all == "no":
        ctl.execute("yhqueue -u `whoami`")

