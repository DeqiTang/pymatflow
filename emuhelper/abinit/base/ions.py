#!/usr/bin/env python
# _*_ coding: utf-8 _*_

class abinit_ions:
    """
    """
    def __init__(self):
        self.params = {
                "ionmov": None,
                "optcell": None,
                "ntime": None,
                "tolmxde": None,
                "tolmxf": None,
                }
    def to_in(self, fout):
        # fout: a file stream for writing
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item ,str(self.params[item])))
        fout.write("\n")
        #
