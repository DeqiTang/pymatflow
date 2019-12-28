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

class cp2k_swarm_global_opt_history:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&HISTORY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END HISTORY\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_swarm_global_opt_minima_crawling_minima_trajectory_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_swarm_global_opt_minima_crawling_minima_trajectory:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_swarm_global_opt_minima_crawling_minima_trajectory_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&MINIMA_TRAJECTORY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END MINIMA_TRAJECTORY\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_swarm_global_opt_minima_crawling:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.minima_trajectory = cp2k_swarm_global_opt_minima_crawling_minima_trajectory()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&MINIMA_CRAWLING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.minima_trajectory.status == True:
            self.minima_trajectory.to_input(fout)
        fout.write("\t\t&END MINIMA_CRAWLING\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "MINIMA_TRAJECTORY":
                slef.minima_trajectory.set_params({item: params[item]})
            else:
                pass


class cp2k_swarm_global_opt_minima_hopping:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&MINIMA_HOPPING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END MINIMA_HOPPING\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_swarm_global_opt_progress_trajectory_each:
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


class cp2k_swarm_global_opt_progress_trajectory:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_swarm_global_opt_progres_trajectory_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&PROGRESS_TRAJECTORY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PROGRESS_TRAJECTORY\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_swarm_global_opt:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.history = cp2k_swarm_global_opt_history()
        self.minima_crawling = cp2k_swarm_global_opt_minima_crawling()
        self.minima_hopping = cp2k_swarm_global_opt_minima_hopping()
        self.progress_trajectory = cp2k_swarm_global_opt_progress_trajectory()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&GLOBAL_OPT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.history.status == True:
            self.history.to_input(fout)
        if self.minima_crawling.status == True:
            self.minima_crawling.to_input(fout)
        if self.minima_hopping.status == True:
            self.minima_hopping.to_input(fout)
        if self.progress_trajectory.status == True:
            self.progress_trajectory.to_input(fout)
        fout.write("\t&END GLOBAL_OPT\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "HISTORY":
                self.history.set_params({item: params[item]})
            elif item.split("-")[1] == "MINIMA_CRAWLING":
                self.minima_crawling.set_params({item: params[item]})
            elif item.split("-")[1] == "MINIMA_HOPPING":
                self.minima_hopping.set_params({item: params[item]})
            elif item.split("-")[1] == "PROGRESS_TRAJECTORY":
                self.progress_trajectory.set_params({item: params[item]})
            else:
                pass


class cp2k_swarm_print_communication_log_each:
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


class cp2k_swarm_print_communication_log:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.each = cp2k_swarm_print_communication_log_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&COMMUNICATION_LOG\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END COMMUNICATION_LOG\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_swarm_print_master_run_info_each:
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


class cp2k_swarm_print_master_run_info:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_swarm_print_master_run_info_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&MASTER_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END MASTER_RUN_INFO\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_swarm_print_worker_run_info_each:
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


class cp2k_swarm_print_worker_run_info:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_swarm_print_worker_run_info_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&WORKER_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END WORKER_RUN_INFO\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_swarm_print:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.communication_log = cp2k_swarm_print_communication_log()
        self.master_run_info = cp2k_swarm_print_master_run_info()
        self.worker_run_info = cp2k_swarm_print_worker_run_info()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.communication_log.status == True:
            self.communication_log.to_input(fout)
        if self.master_run_info.status == True:
            self.master_run_info.to_input(fout)
        if self.worker_run_info.status == True:
            self.worker_run_info.to_input(fout)
        fout.write("\t&END PRINT\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "COMMUNICATION_LOG":
                slef.communication_log.set_params({item: params[item]})
            elif item.split("-")[1] == "MASTER_RUN_INFO":
                self.master_run_info.set_params({item: params[item]})
            elif item.split("-")[1] == "WORKER_RUN_INFO":
                self.worker_run_info.set_params({item: params[item]})
            else:
                pass



class cp2k_swarm:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.global_opt = cp2k_swarm_global_opt()
        self.printout = cp2k_swarm_print()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("&SWARM\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t%s %s\n" % (item, self.params[item]))
        if self.global_opt.status == True:
            self.global_opt.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("&END SWARM\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 1:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[0] == "GLBOAL_OPT":
                self.global_opt.set_params({item: params[item]})
            elif item.split("-")[0] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass


