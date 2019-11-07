#!/usr/bin/env python
# _*_ coding: utf-8 _*_

class cp2k_motion_shell_opt:
    def __init__(self):
        self.params = {
                "MAX_DR": None,
                "MAX_FORCE": None,
                "MAX_ITER": None,
                "RMS_DR": None,
                "RMS_FORCE": None,
                "OPTIMIZER": None,
                "STEP_START_VAL": None,
                }
    def to_motion(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&SHELL_OPT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END SHELL_OPT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]

