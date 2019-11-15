#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ========================
# CP2K / MOTION / CELL_OPT
# ========================
class cp2k_motion_cell_opt:
    def __init__(self):
        self.params = {
                "CONSTRAINT": None,
                "EXTERNAL_PRESSURE": None,
                "KEEP_ANGLES": None,
                "KEEP_SYMMETRY": None,
                "MAX_DR": None,
                "MAX_FORCE": None,
                "MAX_ITER": None,
                "RMS_DR": None,
                "RMS_FORCE": None,
                "OPTIMIZER": None,
                "PRESSURE_TOLLERANCE": None,
                "STEP_START_VAL": None,
                "TYPE": None,  # DIRECT_CELL_OPT, GEO_OPT, MD
                }
        self.default_set()

    def to_motion(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&CELL_OPT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&PRINT\n")
        fout.write("\t\t\t&CELL HIGH\n")
        fout.write("\t\t\t\tFILENAME cell.xyz\n")
        fout.write("\t\t\t&END CELL\n")
        fout.write("\t\t&END PRINT\n")
        fout.write("\t\t&END CELL_OPT\n")
    
    def default_set(self):
        self.params["CONSTRAINT"] = "NONE"
        self.params["KEEP_ANGLES"] = ".FALSE."
        self.params["KEEP_SYMMETRY"] = ".FALSE."
        self.params["MAX_DR"] = 3.0E-3
        self.params["MAX_FORCE"] = 4.5e-4
        self.params["MAX_ITER"] = 200
        self.params["OPTIMIZER"] = "BFGS"
        self.params["PRESSURE_TOLERANCE"] = 1.0E2
        self.params["RMS_DR"] = 1.5e-3
        self.params["RMS_FORCE"] = 3.0e-4
        self.params["TYPE"] = "DIRECT_CELL_OPT"

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]


