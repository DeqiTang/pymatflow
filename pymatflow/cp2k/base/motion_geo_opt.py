#!/usr/bin/env python
# _*_ coding: utf-8 _*_


class cp2k_motion_geo_opt_bfgs_restart_each:
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
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_geo_opt_bfgs_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_geo_opt_bfgs_restart_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END RESTART\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_geo_opt_bfgs:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.restart = cp2k_motion_geo_opt_bfgs_restart()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&BFGS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.restart.status == True:
            self.restart.to_input(fout)
        fout.write("\t\t&END BFGS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "RESTART":
                self.restart.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_geo_opt_cg_line_search_2pnt:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&2PNT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END 2PNT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_geo_opt_cg_line_search_gold:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&GOLD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END GOLD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_geo_opt_cg_line_search:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self._2pnt = cp2k_motion_geo_opt_cg_line_search_2pnt()
        self.gold = cp2k_motion_geo_opt_cg_line_search_gold()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&LINE_SEARCH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self._2pnt.status == True:
            self._2pnt.to_input(fout)
        if self.gold.status == True:
            self.gold.to_input(fout)
        fout.write("\t\t\t&END LINE_SEARCH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "2PNT":
                slef._2pnt.set_params({item: params[item]})
            elif item.split("-")[3] == "GOLD":
                self.gold.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_geo_opt_cg:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.line_serach = cp2k_motion_geo_opt_cg_line_search()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&CG\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.line_search.status == True:
            self.line_search.to_input(fout)
        fout.write("\t\t&END CG\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "LINE_SEARCH":
                self.line_search.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_geo_opt_lbfgs:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&LBFGS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END LBFGS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_geo_opt_print_program_run_info_each:
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
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_geo_opt_print_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.each = cp2k_motion_geo_opt_print_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
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


class cp2k_motion_geo_opt_print:
    def __init__(self):
        self.params = {
                }
        self.status = False
    
        self.program_run_info = cp2k_motion_geo_opt_print_program_run_info()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
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


class cp2k_motion_geo_opt_transition_state_dimer_dimer_vector:
    def __init__(self):
        self.params = {
                }
        self.status = False
    
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&DIMER_VECTOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END DIMER_VECTOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass








class cp2k_motion_geo_opt_transition_state_dimer_rot_opt_bfgs_restart_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_geo_opt_transition_state_dimer_rot_opt_bfgs_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_geo_opt_transition_state_dimer_rot_opt_bfgs_restart_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t&END RESTART\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[6] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_geo_opt_transition_state_dimer_rot_opt_bfgs:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.restart = cp2k_motion_geo_opt_transition_state_dimer_rot_opt_bfgs_restart()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&BFGS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.restart.status == True:
            self.restart.to_input(fout)
        fout.write("\t\t\t\t\t&END BFGS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "RESTART":
                self.restart.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_geo_opt_transition_state_dimer_rot_opt_cg_line_search_2pnt:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&2PNT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t\t&END 2PNT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_geo_opt_transition_state_dimer_rot_opt_cg_line_search_gold:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&GOLD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t\t&END GOLD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_geo_opt_transition_state_dimer_rot_opt_cg_line_search:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self._2pnt = cp2k_motion_geo_opt_transition_state_dimer_rot_opt_cg_line_search_2pnt()
        self.gold = cp2k_motion_geo_opt_transition_state_dimer_rot_opt_cg_line_search_gold()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&LINE_SEARCH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self._2pnt.status == True:
            self._2pnt.to_input(fout)
        if self.gold.status == True:
            self.gold.to_input(fout)
        fout.write("\t\t\t\t\t\t&END LINE_SEARCH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[6] == "2PNT":
                slef._2pnt.set_params({item: params[item]})
            elif item.split("-")[6] == "GOLD":
                self.gold.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_geo_opt_transition_state_dimer_rot_opt_cg:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.line_serach = cp2k_motion_geo_opt_transition_state_dimer_rot_opt_cg_line_search()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&CG\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.line_search.status == True:
            self.line_search.to_input(fout)
        fout.write("\t\t\t\t\t&END CG\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "LINE_SEARCH":
                self.line_search.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_geo_opt_transition_state_dimer_rot_opt_lbfgs:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&LBFGS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END LBFGS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_geo_opt_transition_state_dimer_rot_opt_print_program_run_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_geo_opt_transition_state_dimer_rot_opt_print_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.each = cp2k_motion_geo_opt_transition_state_dimer_rot_opt_print_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_geo_opt_transition_state_dimer_rot_opt_print_rotational_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_geo_opt_transition_state_dimer_rot_opt_print_rotational_info:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.each = cp2k_motion_geo_opt_transition_state_dimer_rot_opt_print_rotational_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&ROTATIONAL_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t&ENDROTATIONAL_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[6] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_geo_opt_transition_state_dimer_rot_opt_print:
    def __init__(self):
        self.params = {
                }
        self.status = False
    
        self.program_run_info = cp2k_motion_geo_opt_transition_state_dimer_rot_opt_print_program_run_info()
        self.rotational_info = cp2k_motion_geo_opt_transition_state_dimer_rot_opt_print_rotational_info()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.rotational_info.status == True:
            self.rotational_info.to_input(fout)
        fout.write("\t\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[5] == "ROTATIONAL_INFO":
                self.rotational_info.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_geo_opt_transition_state_dimer_rot_opt:
    def __init__(self):
        self.params = {
                }
        self.status = False
    
        self.bfgs = cp2k_motion_geo_opt_bfgs()
        self.cg = cp2k_motion_geo_opt_cg()
        self.lbfgs = cp2k_motion_geo_opt_lbfgs()
        self.printout = cp2k_motion_geo_opt_print()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&ROT_OPT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.bfgs.status == True:
            self.bfgs.to_input(fout)
        if self.cg.status == True:
            self.cg.to_input(fout)
        if self.lbfgs.status == True:
            self.lbfgs.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t\t&END ROT_OPT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "BFGS":
                self.bfgs.set_params({item: params[item]})
            elif item.split("-")[4] == "CG":
                self.cg.set_params({item: params[item]})
            elif item.split("-")[4] == "LBFGS":
                self.lbfgs.set_params({item: params[item]})
            elif item.split("-")[4] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_geo_opt_transition_state_dimer:
    def __init__(self):
        self.params = {
                }
        self.status = False
    
        self.dimer_vector = cp2k_motion_geo_opt_transition_state_dimer_dimer_vector()
        self.rot_opt = cp2k_motion_geo_opt_transition_state_dimer_rot_opt()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&DIMER\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.dimer_vector.status == True:
            self.dimer_vector.to_input(fout)
        if self.rot_opt.status == True:
            self.rot_opt.to_input(fout)
        fout.write("\t\t\t&END DIMER\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "DIMER_VECTOR":
                self.dimer_vector.set_params({item: params[item]})
            elif item.split("-")[3] == "ROT_OPT":
                self.rot_opt.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_geo_opt_transition_state:
    def __init__(self):
        self.params = {
                }
        self.status = False
    
        self.dimer = cp2k_motion_geo_opt_transition_state_dimer()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&TRANSITION_STATE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.dimer.status == True:
            self.dimer.to_input(fout)
        fout.write("\t\t&END TRANSITION_STATE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "DIMER":
                self.dimer.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_geo_opt:
    def __init__(self):
        self.params = {
                "MAX_DR": None,
                "MAX_FORCE": None,
                "MAX_ITER": None,
                "RMS_DR": None,
                "RMS_FORCE": None,
                "OPTIMIZER": None, # BFGS(default), CG, LBFGS
                "STEP_START_VAL": None,
                "TYPE": None, # MINIMIZATION(default), TRANSITION_STATE
                }
        self.status = False
        
        self.bfgs = cp2k_motion_geo_opt_bfgs()
        self.cg = cp2k_motion_geo_opt_cg()
        self.lbfgs = cp2k_motion_geo_opt_lbfgs()
        self.printout = cp2k_motion_geo_opt_print()
        self.transition_state = cp2k_motion_geo_opt_transition_state()
        # basic setting
        self.params["MAX_DR"] = 3.0e-3
        self.params["MAX_FORCE"] = 4.5e-4
        self.params["MAX_ITER"] = 200
        self.params["OPTIMIZER"] = "BFGS"
        self.params["RMS_DR"] = 1.5e-3
        self.params["RMS_FORCE"] = 3.0e-4
        self.params["TYPE"] = "MINIMIZATION"


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&GEO_OPT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, str(self.params[item])))
        if self.bfgs.status == True:
            self.bfgs.to_input(fout)
        if self.cg.status == True:
            self.cg.to_input(fout)
        if self.lbfgs.status == True:
            self.lbfgs.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.transition_state.status == True:
            self.transition_state.to_input(fout)
        fout.write("\t&END GEO_OPT\n")
    
         
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "BFGS":
                self.bfgs.set_params({item: params[item]})
            elif item.split("-")[1] == "CG":
                self.cg.set_params({item: params[item]})
            elif item.split("-")[1] == "LBFGS":
                self.lbfgs.set_params({item: params[item]})
            elif item.split("-")[1] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[1] == "TRANSITION_STATE":
                self.transition_state.set_params({item: params[item]})
            else:
                pass
