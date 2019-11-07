#!/usr/bin/env python
# _*_ coding: utf-8 _*_

class cp2k_motion_geo_opt:
    def __init__(self):
        self.params = {
                "MAX_DR": None,
                "MAX_FORCE": None,
                "MAX_ITER": None,
                "RMS_DR": None,
                "RMS_FORCE": None,
                "OPTIMIZER": None, # BFGS(default), CG, LBFGS
                "STEP_START_VAL": None,
                "TYPE": None, # MINIMIZATION(default), TRANSITION_STATE
                }
        self.default_set()

    def to_motion(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&GEO_OPT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t&END GEO_OPT\n")
    
    def default_set(self):
        self.params["MAX_DR"] = 3.0e-3
        self.params["MAX_FORCE"] = 4.5e-4
        self.params["MAX_ITER"] = 200
        self.params["OPTIMIZER"] = "BFGS"
        self.params["RMS_DR"] = 1.5e-3
        self.params["RMS_FORCE"] = 3.0e-4
        self.params["TYPE"] = "MINIMIZATION"
 
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]

