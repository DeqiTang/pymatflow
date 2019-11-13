#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
from emuhelper.qe.static import static_run
from emuhelper.remote.ssh import ssh
from emuhelper.remote.rsync import rsync


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

    # server handle
    if args.auto == 0:
        pass
    elif args.auto == 1:
        mover = rsync()
        mover.get_info(os.path.join(os.path.expanduser("~"), ".emuhelper/server.conf"))
        mover.copy_default(source=os.path.abspath(args.directory))
    elif args.auto == 2:
        mover = rsync()
        mover.get_info(os.path.join(os.path.expanduser("~"), ".emuhelper/server.conf"))
        mover.copy_default(source=os.path.abspath(args.directory))
        ctl = ssh()
        ctl.get_info(os.path.join(os.path.expanduser('~'), ".emuhelper/server.conf"))
        ctl.login()
        ctl.submit(workdir=args.directory, jobfile="relax.in.sub")
