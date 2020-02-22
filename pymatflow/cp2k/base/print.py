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
# CP2K / PRINT
#=============================================

class cp2k_print_distribution_each:
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
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_print_distribution:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_print_distribution_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&distribution\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == true:
            self.each.to_input(fout)
        fout.write("\t\t&end distribution\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "each":
                self.each.set_params({item: params[item]})
            else:
                pass



class cp2k_print_distribution1d_each:
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
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_print_distribution1d:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_print_distribution1d_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&distribution1d\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == true:
            self.each.to_input(fout)
        fout.write("\t\t&end distribution1d\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "each":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_print_distribution2d_each:
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
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_print_distribution2d:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_print_distribution2d_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&distribution2d\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == true:
            self.each.to_input(fout)
        fout.write("\t\t&end distribution2d\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "each":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_print_forces_each:
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
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_print_forces:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_print_forces_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&forces\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == true:
            self.each.to_input(fout)
        fout.write("\t\t&end forces\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "each":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_print_grid_information_each:
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
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_print_grid_information:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_print_grid_information_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&grid_information\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == true:
            self.each.to_input(fout)
        fout.write("\t\t&end grid_information\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "each":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_print_program_run_info_each:
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
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_print_program_run_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_print_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&program_run_info\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == true:
            self.each.to_input(fout)
        fout.write("\t\t&end program_run_info\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "each":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_print_stress_tensor_each:
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
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_print_stress_tensor:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_print_stress_tensor_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&stress_tensor\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == true:
            self.each.to_input(fout)
        fout.write("\t\t&end stress_tensor\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "each":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_print_total_numbers_each:
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
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_print_total_numbers:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_print_total_numbers_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&total_numbers\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == true:
            self.each.to_input(fout)
        fout.write("\t\t&end total_numbers\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "each":
                self.each.set_params({item: params[item]})
            else:
                pass



class cp2k_print:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.distribution = cp2k_print_distribution()
        self.distribution1d = cp2k_print_distribution1d()
        self.distribution2d = cp2k_print_distribution2d()
        self.forces = cp2k_print_forces()
        self.grid_information = cp2k_print_grid_information()
        self.program_run_info = cp2k_print_program_run_info()
        self.stress_tensor = cp2k_print_stress_tensor()
        self.total_numbers = cp2k_print_total_numbers()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.distribution.status == True:
            self.distribution.to_input(fout)
        if self.distribution1d.status == True:
            self.distribution1d.to_input(fout)
        if self.distribution2d.status == True:
            self.distribution2d.to_input(fout)
        if self.forces.status == True:
            self.forces.to_input(fout)
        if self.grid_information.status == True:
            self.grid_information.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.stress_tensor.status == True:
            self.stress_tensor.to_input(fout)
        if self.total_numbers.status == True:
            self.total_numbers.to_input(fout)
        fout.write("\t&END PRINTn")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "DISTRIBUTION":
                self.distribution.set_params({item: params[item]})
            elif item.split("-")[1] == "DISTRIBUTION1D":
                self.distribution1d.set_params({item: params[item]})
            elif item.split("-")[1] == "DISTRIBUTION2D":
                self.distribution2d.set_params({item: params[item]})
            elif item.split("-")[1] == "FORCES":
                self.forces.set_params({item: params[item]})
            elif item.split("-")[1] == "GRID_INFORMATION":
                self.grid_information.set_params({item: params[item]})
            elif item.split("-")[1] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[1] == "STRESS_TENSOR":
                self.stress_tensor.set_params({item: params[item]})
            elif item.split("-")[1] == "TOTAL_NUMBERS":
                self.total_numbers.set_params({item: params[item]})
            else:
                pass
