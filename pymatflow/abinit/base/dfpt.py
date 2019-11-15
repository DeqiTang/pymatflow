#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import numpy as np


class abinit_dfpt:
    """
    see: 
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
                ""
                }
        self.kpoints = kpoints()
                
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
                fout.write("%s %s\n" % (item, str(self.params[item])))
                fout.write("\n")
        fout.write("\n")

    
    def basic_setting(self):
        self.params["rfphon"] = 1
        self.params["rfdir"] = [1, 1, 1]
        self.params["rfelfd"] = 3

    def set_params(self, params):
        for item in params:
            self.params[item] = params[item]

