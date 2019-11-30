#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import numpy as np


class abinit_dfpt:
    """
    Reference:
        https://docs.abinit.org/guide/respfn/
        https://docs.abinit.org/tutorial/rf1/
        https://docs.abinit.org/tutorial/rf2/
        https://docs.abinit.org/tutorial/elastic/
        https://docs.abinit.org/tutorial/ffield/
        https://docs.abinit.org/tutorial/nlo/
        https://docs.abinit.org/topics/DFPT/
        https://docs.abinit.org/variables/dfpt/
    """
    def __init__(self):
        self.params = {
                "rfatpol": None,
                "rfddk": None,
                "rfelfd": None,
                "rfphon": None,
                "rfstrs": None,
                "rfdir": None,
                "rfmagn": None,
                "rfmeth": None,

                "asr": None,
                "dfpt_sciss": None,
                }
        #self.kpoints = kpoints()
                
    def to_in(self, fout):
        # fout: a file stream for writing
        # ------------
        # 检查输入参数
        #self.check_all_params()
        # ---------------
        # 检查输入参数结束
        fout.write("# =====================================\n")
        fout.write("# DFPT related setting\n")
        fout.write("# =====================================\n")
        fout.write("\n")
        for item in self.params:
            if self.params[item] is not None:
                if item == "rfdir":
                    fout.write("rfdir %d %d %d\n" % (self.params[item][0], self.params[item][1], self.params[item][2]))
                elif item == "rfatpol":
                    fout.write("rfatpol %d %d\n" % (self.params["rfatpol"][0], self.params["rfatpol"][1]))
                    pass
                elif item == "qpt":
                    fout.write("qpt %f %f %f\n" % (self.params[item][0], self.params[item][1], self.params[item][2]))
                elif item == "ngqpt":
                    fout.write("ngqpt %d %d %d\n" % (self.params[item][0], self.params[item][1], self.params[item][2]))
                else:
                    fout.write("%s %s\n" % (item, str(self.params[item])))
                fout.write("\n")
        fout.write("\n")

    
    def basic_setting(self):
        self.params["rfdir"] = [1, 1, 1]
        #self.params["nqpt"] = 1
        #self.params["qpt"] = [0, 0, 0]

    def set_params(self, params):
        for item in params:
            self.params[item] = params[item]

