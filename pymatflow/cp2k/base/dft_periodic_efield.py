#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ================================
# ================================
class cp2k_dft_periodic_efield:
    def __init__(self):
        self.section = "FALSE"
        self.params = {
                }

    def to_dft(self, fout):
        fout.write("\t\t&PERIODIC_EFIELD\n")
        for item in self.params:
            if self.params[item] is not None:
                if item == "POLARIZATION":
                    fout.write("\t\t\t%s %.f %.f %.f\n" % (item, self.params[item][0], self.params[item][1], self.params[item][2]))
                else:
                    fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END PERIODIC_EFIELD\n")

