#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

from pymatflow.cp2k.base.atom_print import cp2k_atom_print

"""
Usage:
"""




class cp2k_atom_powell:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&POWELL\n")
        for item in self.params:
            fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END POWELL\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



