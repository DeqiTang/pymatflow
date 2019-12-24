#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os

import argparse
from pymatflow.remote.ssh import ssh


"""
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cmd", type=str, default="pwd",
            help="command line to execute on server")
    args = parser.parse_args()
    # server handle
    ctl = ssh()
    ctl.get_info(os.path.join(os.path.expanduser('~'), ".emuhelper/server.conf"))
    ctl.login()
    ctl.execute(args.cmd)
