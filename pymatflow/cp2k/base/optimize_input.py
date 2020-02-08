#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil


"""
Usage:
"""

class cp2k_optimize_input_force_matching_compare_energies_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_optimize_input_force_matching_compare_energies:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_optimize_input_force_matching_compare_energies_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&COMPARE_ENERGIES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END COMPARE_ENERGIES\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_optimize_input_force_matching_compare_forces_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_optimize_input_force_matching_compare_forces:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_optimize_input_force_matching_compare_forces_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&COMPARE_FORCES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END COMPARE_FORCES\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_optimize_input_force_matching:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.compare_energies = cp2k_optimize_input_force_matching_compare_energies()
        self.compare_forces = cp2k_optimize_input_force_matching_compare_forces()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&FORCE_MATCHING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.compare_eenrgies.status == True:
            self.compare_energies.to_input(fout)
        if self.compare_forces.status == True:
            self.compare_forces.to_input(fout)
        fout.write("\t&END FORCE_MATCHING\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "COMPARE_ENERGIES":
                self.compare_energies.set_params({item: params[item]})
            elif item.split("-")[1] == "COMPARE_FORCES":
                self.compare_forces.set_params({item: params[item]})
            else:
                pass


class cp2k_optimize_input_history_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END EACH\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_optimize_input_history:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_optimize_input_history_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&HISTORY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t&END HISTORY\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_optimize_input_restart_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END EACH\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_optimize_input_restart:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_optimize_input_restart_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t&END RESTART\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_optimize_input_variable:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&VARIABLE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END VARIABLE\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_optimize_input:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.force_matching = cp2k_optimize_input_force_matching()
        self.history = cp2k_optimize_input_history()
        self.restart = cp2k_optimize_input_restart()
        self.variable = cp2k_optimize_input_variable()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("&OPTIMIZE_INPUT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t%s %s\n" % (item, self.params[item]))
        if self.force_matching.status == True:
            self.force_matching.to_input(fout)
        if self.history.status == True:
            self.history.to_input(fout)
        if self.restart.status == True:
            self.restart.to_input(fout)
        if self.variable.status == True:
            self.restart.to_input(fout)
        fout.write("&END OPTIMIZE_INPUT\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 1:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[0] == "FORCE_MATCHING":
                self.force_matching.set_params({item: params[item]})
            elif item.split("-")[0] == "HISTORY":
                self.history.set_params({item: params[item]})
            elif item.split("-")[0] == "RESTART":
                self.restart.set_params({item: params[item]})
            elif item.split("-")[0] == "VARIABLE":
                self.variable.set_params({item: params[item]})
            else:
                pass


