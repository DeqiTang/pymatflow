#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ===================
# ===================
class cp2k_dft_sic:
    def __init__(self):
        self.params = {
                "ORBITAL_SET": None,
                "SIC_METHOD": None,
                "SIC_SCALING_A": None,
                "SIC_SCALING_B": None,
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&SIC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END SIC\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass
