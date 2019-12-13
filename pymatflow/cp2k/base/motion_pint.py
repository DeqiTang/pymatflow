#!/usr/bin/env python
# _*_ coding: utf-8 _*_

class cp2k_motion_pint:
    def __init__(self):
        self.params = {
                "DT": None,
                "FIX_CENTROID_POS": None,
                "HARM_INT": None,
                "ITERATION": None,
                "MAX_STEP": None,
                "NRESPA": None,
                "NUM_STEPS": None,
                "P": None,
                "PROC_PER_REPLICA": None,
                "PROPAGATOR": None,
                "TEMP": None,
                "TRANSFORMATION": None,
                "T_TOL": None,
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&PINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END PINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]


