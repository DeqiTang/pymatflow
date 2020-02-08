#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil


"""
usage:
"""

# ============================================
# CP2K / FORCE_EVAL / EMBED
#=============================================
class cp2k_embed_mapping_force_eval_fragment:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&FRAGMENT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END FRAGMENT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_embed_mapping_force_eval:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.fragment = cp2k_embed_mapping_force_eval_fragment()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&FORCE_EVAL\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.fragment.status == True:
            self.fragment.to_input(fout)
        fout.write("\t\t\t&END FORCE_EVAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "FRAGMENT":
                self.fragment.set_params({item: params[item]})
            else:
                pass

class cp2k_embed_mapping_force_eval_embed_fragment:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&FRAGMENT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END FRAGMENT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_embed_mapping_force_eval_embed:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.fragment = cp2k_embed_mapping_force_eval_embed_fragment()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&FORCE_EVAL_EMBED\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.fragment.status == True:
            self.fragment.to_input(fout)
        fout.write("\t\t\t&END FORCE_EVAL_EMBED\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "FRAGMENT":
                self.fragment.set_params({item: params[item]})
            else:
                pass


class cp2k_embed_mapping:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.force_eval = cp2k_embed_mapping_force_eval()
        self.force_eval_embed = cp2k_embed_mapping_force_eval_embed()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&MAPPING\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.force_eval.status == True:
            self.force_eval.to_input(fout)
        if self.force_eval_embed.status == True:
            self.force_eval_embed.to_input(fout)
        fout.write("\t\t&END MAPPING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "FORCE_EVAL":
                self.force_eval.set_params({item: params[item]})
            elif item.split("-")[2] == "FORCE_EVAL_EMBED":
                self.force_eval_embed.set_params({item: params[item]})
            else:
                pass

class cp2k_embed_print_program_run_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


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


class cp2k_embed_print_program_run_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_embed_print_program_run_info_each()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_embed_print:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.program_run_info = cp2k_embed_print_program_run_info()
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
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass



class cp2k_embed:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.mapping = cp2k_embed_mapping()
        self.printout = cp2k_embed_print()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&EMBED\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.mapping.status == True:
            self.mappint.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t&END DFT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "MAPPING":
                self.mappint.set_params({item: params[item]})
            elif item.split("-")[1] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass
