#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ====================
# CP2K / MOTION / BAND
# ====================
class cp2k_motion_band:
    def __init__(self):
        self.params = {
                "ALIGN_FRAMES": None,
                "BAND_TYPE": None, # CI-NEB, IT-NEB, SM
                "K_SPRING": None,
                "NPROC_REP": None,
                "NUMBER_OF_REPLICA": None,
                "POT_TYPE": None,
                "PROC_DIST_TYPE": None,
                "ROTATE_FRAMES": None,
                "USE_COLVARS": None,
                }
    def to_motion(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&BAND\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t&END BAND\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]


