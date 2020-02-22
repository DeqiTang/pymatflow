#!/usr/bin/env python
# _*_ coding: utf-8 _*_

class cp2k_motion_flexible_partitioning_control_each:
    def __init__(self):
        self.params = {
                }
        self.status = False
        

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END EACH\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_flexible_partitioning_control:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.each = cp2k_motion_flexible_partitioning_control_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&CONTROL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t&END CONTROL\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_flexible_partitioning_weights_each:
    def __init__(self):
        self.params = {
                }
        self.status = False
        

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END EACH\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_flexible_partitioning_weights:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.each = cp2k_motion_flexible_partitioning_weights_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&WEIGHTS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END WEIGHTS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_flexible_partitioning:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.control = cp2k_motion_flexible_partitioning_control()
        self.weights = cp2k_motion_flexible_partitioning_weights()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&FLEXIBLE_PARTITIONING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, str(self.params[item])))
        if self.control.status == True:
            self.control.to_input(fout)
        if self.weights.status == True:
            self.weights.to_input(fout)
        fout.write("\t&END FLEXIBLE_PARTITIONING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "CONTROL":
                self.control.set_params({item: params[item]})
            elif item.split("-")[1] == "WEIGHTS":
                self.weights.set_params({item: params[item]})
            else:
                pass
