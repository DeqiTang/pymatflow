#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ============================
# CP2K / FORCE_EVAL / DFT / QS
# ============================
class cp2k_dft_qs:
    def __init__(self):
        self.params = {
                "METHOD": "GPW",
                "EPS_DEFAULT": "1.0E-10",
                "FORCE_PAW": None,
                }
    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&QS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END QS\n")

    def set_params(self, params):
        for item in params:
            if len(item.spit("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

