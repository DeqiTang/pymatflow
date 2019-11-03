#!/usr/bin/env python
# _*_ coding: utf-8 _*_


class abinit_properties:
    """
    """
    def __init__(self):
        self.params = {
                }

    def to_in(self, fout):
        # fout: a file stream for writing
        fout.write("# ======================================\n")
        fout.write("# properties calculation related setting\n")
        fout.write("# ======================================\n")
        fout.write("\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item, str(self.params[item])))
                fout.write("\n")
        fout.write("\n")
        fout.write("\n")
        #

    def get_option(self, option=[]):
        """
        1: dos
        2: bands
        3: elf
        4:
        """
        if len(option) == 0:
            return
        if 1 in option:
            self.do_dos()
        if 2 in option:
            self.do_bands()
        if 3 in option:
            self.do_elf()
        if 4 in option:
            pass
    
    def do_dos(self):
        self.params["prtdos"] = 1
    
    def do_bands(self):
        self.params["nband"] = 8
        self.params["enunit"] = 1

    def do_elf(self):
        self.params["prtelf"] = 1

    def berry_phase(self):
        self.params["berryopt"] = -1
        self.params["nberry"] = 8
        self.params["dberry"] = "1 1 1 1"
        self.params["rfdir"] = "1 1 1"
         
    def piezoelectric(self):
        self.params["piezoflag"] = 3
