#!/usr/bin/env python

import os
import socket
import configparser
from ssh2.session import Session


class Ssh:
    """
    format of the config file(eg. ~/.emuhelper/server.conf):
    [server]
    ip = xx.xx.xx.xx
    user = xxx
    serverdir = /xxx/xxx/xxx/
    password = xxx
    """
    def __init__(self):
        self.ip = None
        self.user = None
        self.password = None

    def get_info(self, conf):
        config = configparser.ConfigParser()
        config.read(conf)
        self.ip = config.get("server", "ip")
        self.user = config.get("server", "user")
        self.password = config.get("server", "password")
        self.serverdir = config.get("server", "serverdir")

    def login(self):
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.connect((self.ip, 22))

        self.session = Session()
        self.session.handshake(self.sock)
        self.session.userauth_password(self.user, self.password)

    def submit(self, workdir, jobfile, server="pbs"):
        """
        server: pbs or llhpc
        """
        channel = self.session.open_session()
        if server == "llhpc":
            channel.execute("cd %s; yhbatch %s" % (os.path.join(self.serverdir, workdir), jobfile))
        elif server == "pbs":
            channel.execute("cd %s; qsub %s" % (os.path.join(self.serverdir, workdir), jobfile))
        print("\n\n")
        print("=========================================\n")
        print("     information from remote server\n")
        print("=========================================\n")
        all_data = b''
        size, data = channel.read()
        while size > 0:
            #print(data.decode())
            all_data = all_data + data
            size, data = channel.read()
        channel.close()
        print(all_data.decode())
        print("Exit status: %s" % channel.get_exit_status())

    def execute(self, cmd):
        channel = self.session.open_session()
        channel.execute(cmd)
        print("\n\n")
        print("=========================================\n")
        print("     information from remote server\n")
        print("=========================================\n")
        all_data = b''
        size, data = channel.read()
        while size > 0:
            #print(data.decode())
            all_data = all_data + data
            size, data = channel.read()
        channel.close()
        print(all_data.decode())
        print("Exit status: %s" % channel.get_exit_status())

    def execute_silent(self, cmd):
        channel = self.session.open_session()
        channel.execute(cmd)
        all_data = b''
        size, data = channel.read()
        while size > 0:
            all_data = all_data + data
            size, data = channel.read()
        channel.close()
        return all_data
