#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

import pymatgen as mg

"""
usage:
"""

class qe_cell:
    """

    """
    def __init__(self):
        self.params = {
                }

    def to_in(self, fout):
        # fout: a file stream for writing
        fout.write("&cell\n")
        for item in self.params:
            if self.params[item] is not None:
                if type(self.params[item]) == str:
                    fout.write("%s = '%s'\n" % (item, str(self.params[item])))
                else:
                    fout.write("%s = %s\n" % (item, str(self.params[item])))
        fout.write("/\n")
        fout.write("\n")

    def set_params(self, params):
        """
        params: a dict storing the parameters and values
        """
        for item in params:
            self.params[item] = params[item]
