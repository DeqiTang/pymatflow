#!/usr/bin/env python
# _*_ coding: utf-8 _*_

class cp2k_motion_print_trajectory:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t&TRAJECTORY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END TRAJECTORY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_print_restart_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_print_restart:
    def __init__(self):
        self.params = {}
        self.status = False
        self.each = cp2k_motion_print_restart_each()

    def to_input(self, fout):
        fout.write("\t\t&RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        self.each.to_input(fout)
        fout.write("\t\t&END RESTART\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})

class cp2k_motion_print_restart_history_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_print_restart_history:
    def __init__(self):
        self.params = {}
        self.status = False
        self.each = cp2k_motion_print_restart_history_each()

    def to_input(self, fout):
        fout.write("\t\t&RESTART_HISTORY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        self.each.to_input(fout)
        fout.write("\t\t&END RESTART_HISTORY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})

class cp2k_motion_print:
    def __init__(self):
        self.status = False
        self.trajectory = cp2k_motion_print_trajectory()
        self.restart = cp2k_motion_print_restart()
        self.restart_history = cp2k_motion_print_restart_history()
    
    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&PRINT\n")

        self.trajectory.to_input(fout)
        self.restart.to_input(fout)
        self.restart_history.to_input(fout)
        
        fout.write("\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "TRAJECTORY":
                self.trajectory.set_params({item: params[item]})
            elif item.split("-")[1] == "RESTART":
                self.restart.set_params({item: params[item]})
            elif item.split("-")[1] == "RESTART_HISTORY":
                self.restart_history.set_params({item: params[item]})

