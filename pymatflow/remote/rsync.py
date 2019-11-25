#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import configparser

class rsync:
    def __init__(self):
        pass
    def get_info(self, conf):
        config = configparser.ConfigParser()
        config.read(conf)
        self.user = config.get("server", "user")
        self.ip = config.get("server", "ip")
        self.password = config.get("server", "password")
        self.serverdir = config.get("server", "serverdir")

    def copy(self, source, target):
        print("=========================================\n")
        print("command line to execute:\n")
        print("  rsync -av --progress %s %s\n" % (source, target))
        print("......\n")
        os.system("rsync -av %s %s" % (source, target))

    def copy_default(self, source):
        self.copy(source=source, target=self.user+"@"+self.ip+":"+self.serverdir)
