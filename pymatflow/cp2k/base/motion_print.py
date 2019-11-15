#!/usr/bin/env python
# _*_ coding: utf-8 _*_

class cp2k_motion_print_trajectory:
    def __init__(self):
        self.params = {}

    def to_print(self, fout):
        pass

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_print:
    def __init__(self):
        self.trajectory = cp2k_motion_print_trajectory()
    
    def to_motion(self, fout):
        # fout: a file stream for writing
        fout.write("\t&PRINT\n")
        fout.write("\t\t&TRAJECTORY\n")
        for item in self.trajectory.params:
            fout.write("\t\t\t%s %s\n" % (item, self.trajectory.params[item]))
        fout.write("\t\t&END TRAJECTORY\n")
        fout.write("\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "TRAJECTORY":
                self.trajectory.set_params({item: params[item]})

