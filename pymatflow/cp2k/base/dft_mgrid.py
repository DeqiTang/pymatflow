#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ===============================
# CP2K / FORCE_EVAL / DFT / MGRID
# ===============================
class cp2k_dft_mgrid_interpolator_conv_info_each:
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
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&end EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_mgrid_interpolator_conv_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_mgrid_interpolator_conv_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&CONV_INFO\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&end CONV_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_mgrid_interpolator:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.conv_info = cp2k_dft_mgrid_interpolator_conv_info()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&interpolator\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.conv_info.status == True:
            self.conv_info.to_input(fout)
        fout.write("\t\t\t&end interpolator\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "CONV_INFO":
                self.conv_info.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_mgrid_rs_grid:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RS_GRID\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END RS_GRID\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_mgrid:
    def __init__(self):
        self.params = {
                "CUTOFF": None,
                "REL_CUTOFF": None,
                "NGRIDS": 4,
                "COMMENSURATE": None,
                "MULTIGRID_CUTOFF": None,
                "MULTIGRID_SET": None,
                "PROGRESSION_FACTOR": None,
                "SKIP_LOAD_BALANCE_DISTRIBUTED": None,
                }
        self.status = False

        self.interpolator = cp2k_dft_mgrid_interpolator()
        self.rs_grid = cp2k_dft_mgrid_rs_grid()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&MGRID\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.interpolator.status == True:
            self.interpolator.to_input(fout)
        if self.rs_grid.status == True:
            self.rs_grid.to_input(fout)
        fout.write("\t\t&END MGRID\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "INTERPOLATOR":
                self.interpolator.set_params({item: params[item]})
            elif item.split("-")[2] == "RS_GRID":
                self.rs_grid.set_params({item: params[item]})
            else:
                pass

