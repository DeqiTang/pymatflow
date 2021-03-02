#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import configparser

class Rsync:
    def __init__(self):
        pass
    def get_info(self, conf):
        config = configparser.ConfigParser()
        config.read(conf)
        self.user = config.get("server", "user")
        self.ip = config.get("server", "ip")
        self.password = config.get("server", "password")
        self.serverdir = config.get("server", "serverdir")


    def copy(self, source, target, exclude=[]):
        """
        through exlude we can specify files or dirs to exclude from syncing.
        like exclude=['/a', '/b'], it will exclude 'a' and 'b' in root directory, namely source.
        """
        print("=========================================\n")
        print("command line to execute:\n")
        print("  rsync -av --progress %s %s" % (source, target))
        if len(exclude) != 0:
            for i in range(len(exclude)):
                print(" --exclude=%s" % exclude[i])
            print("\n")
        print("......\n")
        cmd = "rsync -av --progress %s %s" % (source, target)
        if len(exclude) != 0:
            for i in range(len(exclude)):
                cmd = cmd + " --exclude=%s" % exclude[i]
        os.system(cmd)

    def copy_default(self, source):
        self.copy(source=source, target=self.user+"@"+self.ip+":"+self.serverdir)
