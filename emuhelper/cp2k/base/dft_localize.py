#!/usr/bin/env python
# _*_ coding: utf-8 _*_


# =================================
# CP2K / FORCE_EVAL /DFT / LOCALIZE
# =================================
class cp2k_dft_localize:
    """
    About:
        to calculate IR spectra from md running
        it is necessary to have dipole information
        for the molecules available in the simulated
        trajectory.
    """
    def __init__(self):
        self.params = {
                "METHOD": None,
                "MAX_ITER": None,
                }
        self.status = False

    def to_dft(self, fout):
        # fout: a file stream for writing
        if self.status == False:
            return
        fout.write("\t\t&LOCALIZE %s\n" % ".TRUE.")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&PRINT\n")
        fout.write("\t\t\t\t&WANNIER_CENTERS\n")
        fout.write("\t\t\t\t\tIONS+CENTERS\n")
        fout.write("\t\t\t\t\tFILENAME xxx_wannier.xyz\n")
        fout.write("\t\t\t\t&END WANNIER_CENTERS\n")
        fout.write("\t\t\t&END PRINT\n")
        fout.write("\t\t&END LOCALIZE\n")


