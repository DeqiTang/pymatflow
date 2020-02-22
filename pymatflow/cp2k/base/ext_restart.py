#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil



"""
Usage:
"""

class cp2k_ext_restart:
    """
    """
    def __init__(self, project_name="ab-initio"):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("&EXT_RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t%s %s\n" % (item, self.params[item]))
        fout.write("&END EXT_RESTART\n")
        fout.write("\n")
    
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 1:
                self.params[item.split("-")[-1]] = params[item]
