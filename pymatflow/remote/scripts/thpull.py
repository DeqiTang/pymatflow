#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import sys

import argparse
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync


"""
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", type=str,
            help="directory to pull from origin to local")

    parser.add_argument("--tmp", type=str, default="no",
            choices=["yes", "no"],
            help="whether to pull xx/tmp directory back, default no")

    args = parser.parse_args()
    # server handle

    mover = rsync()
    mover.get_info(os.path.join(os.path.expanduser("~"), ".emuhelper/server.conf"))
    
    if args.tmp == "no":
        ctl = ssh()
        ctl.get_info(os.path.join(os.path.expanduser('~'), ".emuhelper/server.conf"))
        ctl.login()
        dirs = ctl.execute_silent("ls")[1].decode().split()
        if args.directory not in dirs:
            print("====================================\n")
            print("           WARNING !!!\n")
            print("====================================\n")
            print("directory not found in remote server\n")
            print("make sure you type in the right one\n")
            sys.exit()
        ctl.execute_silent("mv %s/tmp ./%s.tmp" % (args.directory, args.directory))
        mover.copy(source=os.path.join(mover.user+"@"+mover.ip+":"+mover.serverdir, args.directory)+"/", target=args.directory)
        ctl.execute_silent("mv %s.tmp %s/tmp" % (args.directory, args.directory))
    elif args.tmp == "yes":
        mover.copy(source=os.path.join(mover.user+"@"+mover.ip+":"+mover.serverdir, args.directory)+"/", target=args.directory)
