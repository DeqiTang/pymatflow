#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ==================================
# ==================================
class cp2k_dft_real_time_propagation_print_current_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_real_time_propagation_print_current:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_real_time_propagation_print_current_each()
        
        # basic setting

    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&CURRENT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END CURRENT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_real_time_propagation_print_program_run_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False
    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_real_time_propagation_print_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_real_time_propagation_print_program_run_info_each()
        
        # basic setting

    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_real_time_propagation_print_restart_each:
    def __init__(self):
        self.params = {
                }
        self.status = False
    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_real_time_propagation_print_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_real_time_propagation_print_restart_each()
        
        # basic setting

    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END RESTART\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item] 
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_real_time_propagation_print_restart_history_each:
    def __init__(self):
        self.params = {
                }
        self.status = False
    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_real_time_propagation_print_restart_history:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_real_time_propagation_print_restart_history_each()
        
        # basic setting

    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&RESTART_HISTORY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END RESTART_HISTORY")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_real_time_propagation_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.current = cp2k_dft_real_time_propagation_print_current()
        self.program_run_info = cp2k_dft_real_time_propagation_print_program_run_info()
        self.restart = cp2k_dft_real_time_propagation_print_restart()
        self.restart_history = cp2k_dft_real_time_propagation_print_restart_history()

        # basic settign
    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.current.status == True:
            self.current.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.restart.status == True:
            self.restart.to_input(fout)
        if self.restart_history.status == True:
            self.restart_history.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "CURRENT":
                self.current.set_params({item: params[item]})
            elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[3] == "RESTART":
                self.restart.set_params({item: params[item]})
            elif item.split("-")[3] == "RESTART_HISTORY":
                self.restart_history.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_real_time_propagation:
    def __init__(self):
        self.params = {
                "ACCURACY_REFINEMENT": None,
                "APPLY_DELTA_PULSE": None,
                "ASPC_ORDER": None,
                "PERIODIC": None,
                "PROPAGATOR": None,
                "MAX_EXP": None,
                "MAX_ITER": None,
                "EPS_ITER": None,
                }
        self.status = False

        self.printout = cp2k_dft_real_time_propagation_print()
        
        # basic setting
    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&REAL_TIME_PROPAGATION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t&END REAL_TIME_PROPAGATION\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass
