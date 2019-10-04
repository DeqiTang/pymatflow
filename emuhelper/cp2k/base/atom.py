#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

"""
Usage:
"""

class cp2k_atom:
    """

    """
    def __init__(self):
        self.params = {
                }

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("&ATOM\n")
        for item in self.params:
            fout.write("\t%s %s\n" % (item, self.params[item]))
        fout.write("&END ATOM\n")
        fout.write("\n")

