#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# =================================
# CP2K / FORCE_EVAL / DFT / POISSON
# =================================
class cp2k_dft_poisson:
    def __init__(self):
        self.params = {
                "PERIODIC": "XYZ",
                "POISSON_SOLVER": "PERIODIC",
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&POISSON\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END POISSON\n")


