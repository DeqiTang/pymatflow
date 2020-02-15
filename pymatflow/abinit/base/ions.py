#!/usr/bin/env python
# _*_ coding: utf-8 _*_

class abinit_ions:
    """
    """
    def __init__(self):
        self.params = {}
        self.incharge = [
            "ionmov", "optcell", "ntime", "tolmxde", "tolmxf",
            ]
    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("# ============================\n")
        fout.write("# ions moving related settting\n")
        fout.write("# ============================\n")
        fout.write("\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item ,str(self.params[item])))
                fout.write("\n")
        fout.write("\n\n")
        #

    def basic_setting(self, mode="opt"):
        """
        mode: opt or md
        """
        if mode == "opt":
            self.params["ionmov"] = 3
            self.params["optcell"] = 0
            self.params["ntime"] = 100
            #self.params["nctime"] = 1 # for netcdf writing
            self.params["tolmxf"] = 0 #5.0e-4 # Ha/Bohr
            self.params["tolmxde"] = "0.0001 eV"
        elif mode == "md":
            self.params["ionmov"] = 8 #11
            self.params["dtion"] = 100
            self.params["ntime"] = 1000
            self.params["nctime"] = 1
            self.params["mdtemp(1)"] = 300
            self.params["mdtemp(2)"] = 300
            self.params["tolmxf"] = 5.0e-4 # Ha/Bohr

    def set_params(self, params):
        for item in params:
            self.params[item] = params[item]
