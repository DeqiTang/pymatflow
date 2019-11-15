#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ===============================
# CP2K / FORCE_EVAL / DFT / MGRID
# ===============================
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

    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&MGRID\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END MGRID\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

