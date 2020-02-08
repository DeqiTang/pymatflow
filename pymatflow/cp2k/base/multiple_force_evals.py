#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil


"""
Usage:
"""


class cp2k_multiple_force_evals:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("&MULTIPLE_FORCE_EVALS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t%s %s\n" % (item, self.params[item]))
        fout.write("&END MULTIPLE_FORCE_EVALS\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 1:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


