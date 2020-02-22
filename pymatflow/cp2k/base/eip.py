#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil


# ============================================
# CP2K / EIP
#=============================================

class cp2k_eip_print_coord_avg_each:
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

class cp2k_eip_print_coord_avg:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_eip_print_coord_avg_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&COORD_AVG\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END COORD_AVG\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_eip_print_coord_var_each:
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


class cp2k_eip_print_coord_var:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_eip_print_coord_var_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&COORD_VAR\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END COORD_VAR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_eip_print_count_each:
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


class cp2k_eip_print_count:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_eip_print_count_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&COUNT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END COUNT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_eip_print_energies_each:
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


class cp2k_eip_print_energies:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_eip_print_energies_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&ENERGIES\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END ENERGIES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_eip_print_energies_var_each:
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


class cp2k_eip_print_energies_var:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_eip_print_energies_var_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&ENERGIES_VAR\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END ENERGIES_VAR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_eip_print_forces_each:
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



class cp2k_eip_print_forces:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_eip_print_forces_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&FORCES\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END FORCES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_eip_print:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.coord_avg = cp2k_eip_print_coord_avg()
        self.coord_var = cp2k_eip_print_coord_var()
        self.count = cp2k_eip_print_count()
        self.energies = cp2k_eip_print_energies()
        self.energies_var = cp2k_eip_print_energies_var()
        self.forces = cp2k_eip_print_forces()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.coord_avg.status == True:
            self.coord_avg.to_input(fout)
        if self.coord_var.status == True:
            self.coord_var.to_input(fout)
        if self.count.status == True:
            self.count.to_input(fout)
        if self.energies.status == True:
            self.energies.to_input(fout)
        if self.energies_var.status == True:
            self.energies_var.to_input(fout)
        if self.forces.status == True:
            self.forces.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "COORD_AVG":
                self.coord_avg.set_params({item: params[item]})
            elif item.split("-")[2] == "COORD_VAR":
                self.coord_var.set_params({item: params[item]})
            elif item.split("-")[2] == "COUNT":
                self.count.set_params({item: params[item]})
            elif item.split("-")[2] == "ENERGIES":
                self.energies.set_params({item: params[item]})
            elif item.split("-")[2] == "ENERGIES_VAR":
                self.energies_var.set_params({item: params[item]})
            elif item.split("-")[2] == "FORCES":
                self.forces.set_params({item: params[item]})
            else:
                pass


class cp2k_eip:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.printout = cp2k_eip_print()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&EIP\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t&END EIP\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass
