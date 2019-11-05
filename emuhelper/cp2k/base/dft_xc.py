#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# =============================
# CP2K / FORCE_EVAL / DFT / XC
# =============================
class cp2k_dft_xc_xc_functional:
    def __init__(self):
        self.section = "PBE"
        self.params = {
                }
    def to_xc(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&XC_FUNCTIONAL %s\n" % self.section)
        fout.write("\t\t\t&END XC_FUNCTIONAL\n")

    def set_params(self, params):
        """
        set_params for xc_functional is different from many other
        set_params, as it deal with the key 'XC_FUNCTIONAL' only
        """
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc:
    def __init__(self):
        self.params = {
                "DENSITY_CUTOFF": None,
                "DENSITY_SMOOTH_CUTOFF_RANGE": None,
                "FUNCTIONAL_ROUTINE": None,
                "GRADIENT_CUTOFF": None,
                "TAU_CUTOFF": None,
                }
        self.xc_functional = cp2k_dft_xc_xc_functional()

    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&XC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s" % (item, str(self.params[item])))
        self.xc_functional.to_xc(fout)
        fout.write("\t\t&END XC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                if item.split("-")[-1] == "XC_FUNCTIONAL":
                    self.xc_functional.section = params[item]
                else:
                    self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "XC_FUNCTIONAL" and len(item.split("-")) > 3:
                self.xc_functional.set_params({item: params[item]})

