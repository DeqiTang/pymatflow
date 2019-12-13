#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# =======================================
# CP2K / FORCE_EVAL / DFT / LOW_SPIN_ROKS
# =======================================
class cp2k_dft_low_spin_roks:
    def __init__(self):
        self.params = {
                "ENERGY_SCALING": None,
                "SPIN_CONFIGURATION": None
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&LOW_SPIN_ROKS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END LOW_SPIN_ROKS\n")


