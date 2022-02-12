#!/usr/bin/env python
# _*_ coding: utf-8 _*_



class cp2k_motion_tmc_move_type:
    def __init__(self):
        self.params = {
                }

        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&MOVE_TYPE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END MOVE_TYPE\n")
    
         
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_tmc_nmc_moves_move_type:
    def __init__(self):
        self.params = {
                }

        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&MOVE_TYPE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END MOVE_TYPE\n")
    
         
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_tmc_nmc_moves:
    def __init__(self):
        self.params = {
                }

        self.status = False

        self.mvoe_type = cp2k_motion_tmc_nmc_moves_move_type()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&NMC_MOVES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.mvoe_type.status == True:
            self.move_type.to_input(fout)
        fout.write("\t\t&END NMC_MOVES\n")
    
         
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "MOVE_TYPE":
                self.move_type.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_tmc_tmc_analysis_charge:
    def __init__(self):
        self.params = {
                }

        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CHARGE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END CHARGE\n")
    
         
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_tmc_tmc_analysis:
    def __init__(self):
        self.params = {
                }

        self.status = False

        self.charge = cp2k_motion_tmc_tmc_analysis_charge()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&TMC_ANALYSIS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.charge.status == True:
            self.charge.to_input(fout)
        fout.write("\t\t&END TMC_ANALYSIS\n")
    
         
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "CHARGE":
                self.charge.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_tmc_tmc_analysis_files_charge:
    def __init__(self):
        self.params = {
                }

        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CHARGE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END CHARGE\n")
    
         
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_tmc_tmc_analysis_files:
    def __init__(self):
        self.params = {
                }

        self.status = False

        self.charge = cp2k_motion_tmc_tmc_analysis_files_charge()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&TMC_ANALYSIS_FILES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.charge.status == True:
            self.charge.to_input(fout)
        fout.write("\t\t&END TMC_ANALYSIS_FILES\n")
    
         
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "CHARGE":
                self.charge.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_tmc:
    def __init__(self):
        self.params = {
                }

        self.status = False

        self.move_type = cp2k_motion_tmc_move_type()
        self.nmc_moves = cp2k_motion_tmc_nmc_moves()
        self.tmc_analysis = cp2k_motion_tmc_tmc_analysis()
        self.tmc_analysis_files = cp2k_motion_tmc_tmc_analysis_files()
        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&TMC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, str(self.params[item])))
        if self.move_type.status == True:
            self.move_type.to_input(fout)
        if self.nmc_moves.status == True:
            self.nmc_moves.to_input(fout)
        if self.tmc_analysis.status == True:
            self.tmc_analysis.to_input(fout)
        if self.tmc_analysis_files.status == True:
            self.tmc_analysis_files.to_input(fout)
        fout.write("\t&END TMC\n")
    
         
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "MOVE_TYPE":
                self.move_type.set_params({item: params[item]})
            elif item.split("-")[1] == "NMC_MOVES":
                self.nmc_moves.set_params({item: params[item]})
            elif item.split("-")[1] == "TMC_ANALYSIS":
                self.tmc_analysis.set_params({item: params[item]})
            elif item.split("-")[1] == "TMC_ANALYSIS_FILES":
                self.tmc_analysis_files.set_params({item: params[item]})
            else:
                pass
