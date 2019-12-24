#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ================================
# ================================
class cp2k_dft_density_fitting_program_run_info_each:
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
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_density_fitting_program_run_info:
    def __init__(self):
        self.params = {

                }
        self.status = False

        self.each = cp2k_dft_density_fitting_program_run_info_each()

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_density_fitting:
    def __init__(self):
        self.params = {

                }
        self.status = False
            
        self.program_run_info = cp2k_dft_density_fitting_program_run_info()
        
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&DENSITY_FITTING\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t&END DENSITY_FITTING\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass
