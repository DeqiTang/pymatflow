#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


"""
Usage:
"""


class cp2k_farming_job:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&JOB\n")
        for item in self.params:
            fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END JOB\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_farming_program_run_info_each:
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
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END EACH\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_farming_program_run_info:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.each = cp2k_farming_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t&END PROGRAM_RUN_INFO\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_farming_restart_each:
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
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END EACH\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_farming_restart:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.each = cp2k_farming_restart_each()
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


class cp2k_farming:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.job = cp2k_farming_job()
        self.program_run_info = cp2k_farming_program_run_info()
        self.restart = cp2k_farming_restart()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("&FARMING\n")
        for item in self.params:
            fout.write("\t%s %s\n" % (item, self.params[item]))
        if self.job.status == True:
            self.job.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.restart.status == True:
            self.restart.to_input(fout)
        fout.write("&END FARMING\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 1:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[0] == "JOB":
                self.job.set_params([item: params[item]])
            elif item.split("-")[0] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[0] == "RESTART":
                self.restart.set_params({item: params[item]})
            else:
                pass


