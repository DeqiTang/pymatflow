#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os

import argparse
from pymatflow.remote.ssh import ssh


"""
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--id", type=str,
            help="JOBID of the task to cancel")
    args = parser.parse_args()
    # server handle
    ctl = ssh()
    ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_yh.conf"))
    ctl.login()
    ctl.execute("yhcancel %s" % args.id)

