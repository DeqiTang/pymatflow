#!/usr/bin/env python
# _*_ coding: utf-8 _*_


class abinit_properties:
    """
    """
    def __init__(self):
        self.params = {
                }

    def to_in(self, fout):
        # fout: a file stream for writing
        fout.write("# properties\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item, str(self.params[item])))
        fout.write("\n")
        #

    def berry_phase(self):
        self.params["berryopt"] = -1
        self.params["nberry"] = 8
        self.params["dberry"] = "1 1 1 1"
        self.params["rfdir"] = "1 1 1"
         
    def piezoelectric(self):
        self.params["piezoflag"] = 3
