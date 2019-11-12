#!/usr/bin/env python

import os
import socket
import configparser
from ssh2.session import Session


class ssh:
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

    def submit(self, workdir, jobfile):
        channel = self.session.open_session()
        channel.execute("cd %s; yhbatch -p free %s" % (os.path.join(self.serverdir, workdir), jobfile))
        print("\n\n")
        print("=========================================\n")
        print("     information from remote server\n")
        print("=========================================\n")
        size, data = channel.read()
        while size > 0:
            print(data.decode())
            size, data = channel.read()
        channel.close()
        print("Exit status: %s" % channel.get_exit_status())

    def execute(self, cmd):
        channel = self.session.open_session()
        channel.execute(cmd)
        print("\n\n")
        print("=========================================\n")
        print("     information from remote server\n")
        print("=========================================\n")
        size, data = channel.read()
        while size > 0:
            print(data.decode())
            size, data = channel.read()
        channel.close()
        print("Exit status: %s" % channel.get_exit_status())

    def execute_silent(self, cmd):
        channel = self.session.open_session()
        channel.execute(cmd)
        #size, data = channel.read()
        #while size > 0:
        #    print(data.decode())
        #    size, data = channel.read()
        out = channel.read()
        # out is a tuple of length 2
        # out[0] is the size of the bytes stored inf out[1]
        channel.close()
        return out
