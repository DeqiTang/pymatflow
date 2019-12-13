#!/usr/bin/env python
# _*_ coding: utf-8 _*_

class cp2k_motion_constraint:
    def __init__(self):
        self.params = {
                "CONSTRAINT_INIT": None,
                "ROLL_TOLERANCE": None,
                "SHAKE_TOLERANCE": None,
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&CONSTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END CONSTRAINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]


