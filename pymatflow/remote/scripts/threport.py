#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os

from pymatflow.remote.ssh import ssh


"""
"""

if __name__ == "__main__":
    # server handle
    ctl = ssh()
    ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_yh.conf"))
    ctl.login()
    ctl.execute("yhreport CLuster UserUtilizationByAccount start=11/1/14 end=now -t hour")

