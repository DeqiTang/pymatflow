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
# CP2K / MIXED
#=============================================


class cp2k_mixed_coupling:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&coupling\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&end coupling\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mixed_generic:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&GENERIC\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&end GENERIC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_mixed_linear:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&LINEAR\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&end LINEAR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mixed_mapping_force_eval_fragment:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&FRAGMENT\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&end FRAGMENT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_mixed_mapping_force_eval:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.fragment = cp2k_mixed_mapping_force_eval_fragment()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&FORCE_EVAL\n")
        for item in self.params:
            if self.params[item] is not none:                    
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


class cp2k_mixed_mapping_force_eval_mixed_fragment:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&FRAGMENT\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&end FRAGMENT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_mixed_mapping_force_eval_mixed:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.fragment = cp2k_mixed_mapping_force_eval_mixed_fragment()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&FORCE_EVAL_MIXED\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.fragment.status == True:
            self.fragment.to_input(fout)
        fout.write("\t\t\t&end FORCE_EVAL_MIXED\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "FRAGMENT":
                self.fragment.set_params({item: params[item]})
            else:
                pass



class cp2k_mixed_mapping:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.force_eval = cp2k_mixed_mapping_force_eval()
        self.force_eval_mixed = cp2k_mixed_mapping_force_eval_mixed()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&MAPPING\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.force_eval.status == True:
            self.force_eval.to_input(fout)
        if self.force_eval_mixed.status == True:
            self.force_eval_mixed.to_input(fout)
        fout.write("\t\t&end MAPPING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "FORCE_EVAL":
                self.force_eval.set_params({item: params[item]})
            elif item.split("-")[2] == "FORCE_EVAL_MIXED":
                self.force_eval_mixed.set_params({item: params[item]})
            else:
                pass



class cp2k_mixed_mixed_cdft_block_diagonalize:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&BLOCK_DIAGONALIZE\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END BLOCK_DIAGONALIZE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mixed_mixed_cdft_print_program_run_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mixed_mixed_cdft_print_program_run_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_mixed_mixed_cdft_print_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_mixed_mixed_cdft_print:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.program_run_info = cp2k_mixed_mixed_cdft_print_program_run_info()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass


class cp2k_mixed_mixed_cdft:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.block_diagonalize = cp2k_mixed_mixed_cdft_block_diagonalize()
        self.printout = cp2k_mixed_mixed_cdft_print()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&MIXED_CDFT\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.block_diagonalize.status == True:
            self.block_diagonalize.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t&end MIXED_CDFT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "BLOCK_DIAGONALIZE":
                self.block_diagonalize.set_params({item: params[item]})
            elif item.split("-")[2] == "PRITN":
                self.printout.set_params({item: params[item]})
            else:
                pass


class cp2k_mixed_print_dipole_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&end EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_mixed_print_dipole:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_mixed_print_dipole_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&end PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_mixed_print_program_run_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&end EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_mixed_print_program_run_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_mixed_print_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass



class cp2k_mixed_print:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.dipole = cp2k_mixed_print_dipole()
        self.program_run_Info = cp2k_mixed_print_program_run_info()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.dipole.status == True:
            self.dipole.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t&end PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "DIPOLE":
                self.dipole.set_params({item: params[item]})
            elif item.split("-")[2] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass


class cp2k_mixed_restraint:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&RESTRAINT\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&end RESTRAINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_mixed:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.coupling = cp2k_mixed_coupling()
        self.generic = cp2k_mixed_generic()
        self.linear = cp2k_mixed_linear()
        self.mapping = cp2k_mixed_mapping()
        self.mixed_cdft  = cp2k_mixed_cdft()
        self.printout = cp2k_mixed_print()
        self.restraint = cp2k_mixed_restraint()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&MIXED\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.coupling.status == True:
            self.coupling.to_input(fout)
        if self.generic.status == True:
            self.generic.to_input(fout)
        if self.linear.status == True:
            self.linear.to_input(fout)
        if self.mapping.status == True:
            self.mapping.to_input(fout)
        if self.mixed_cdft.status == True:
            self.mixed_cdft.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.restraint.status == True:
            self.restraint.to_input(fout)
        fout.write("\t&END MIXED\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "COUPLING":
                self.coupling.set_params({item: params[item]})
            elif item.split("-")[1] == "GENERIC":
                self.generic.set_params({item: params[item]})
            elif item.split("-")[1] == "LINEAR":
                self.linear.set_params({item: params[item]})
            elif item.split("-")[1] == "MAPPING":
                self.mapping.set_params({item: params[item]})
            elif item.split("-")[1] == "MIXED_CDFT":
                self.mixed_cdft.set_params({item: params[item]})
            elif item.split("-")[1] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[1] == "RESTRAITN":
                self.restraint.set_params({item: params[item]})
            else:
                pass
