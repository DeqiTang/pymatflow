#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


"""
usage:
"""

# ============================================
# CP2K / BSSE
#=============================================


class cp2k_bsse_configuration:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&CONFIGURATION\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END CONFIGURATION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_bsse_fragment:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&FRAGMENT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END FRAGMENT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_bsse_fragment_energies:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&FRAGMENT_ENERGIES\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END FRAGMENT_ENERGIES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_bsse_print_program_run_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_bsse_print_program_run_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_bsse_print_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_bsse_print_restart_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_bsse_print_restart:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_bsse_print_restart_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RESTART\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RESTART\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_bsse_print:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.program_run_info = cp2k_bsse_print_program_run_info()
        self.restart = cp2k_bsse_print_restart()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.restart.status == True:
            self.restart.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[2] == "RESTART":
                self.restart.set_params({item: params[item]})
            else:
                pass


class cp2k_bsse:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.configuration = cp2k_bsse_configuration()
        self.fragment = cp2k_bsse_fragment()
        self.fragment_energies = cp2k_bsse_fragment_energies()
        self.printout = cp2k_bsse_print()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&BSSE\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.configuration.status == True:
            self.configuration.to_input(fout)
        if self.fragment.status == True:
            self.fragment.to_input(fout)
        if self.fragment_energies.status == True:
            self.fragment_energies.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t&END BSSE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "CONFIGURATION":
                self.configuration.set_params({item: params[item]})
            elif item.split("-")[1] == "FRAGMENT":
                self.fragment.set_params({item: params[item]})
            elif item.split("-")[1] == "FRAGMENT_ENERGIES":
                self.fragment_energies.set_params({item: params[item]})
            elif item.split("-")[1] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass
