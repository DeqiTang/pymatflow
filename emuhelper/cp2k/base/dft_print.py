#!/usr/bin/env python
# _*_ coding: utf-8 _*_


# ===============================
# CP2K / FORCE_EVAL / DFT / PRINT
# ===============================
class cp2k_dft_print:
    def __init__(self):
        self.params = {
                }
        self.electron_density = False
        self.bands = False
        self.pdos = False
        self.moments = False

    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.electron_density == True:
            fout.write("\t\t\t&E_DENSITY_CUBE\n")
            fout.write("\t\t\t&END E_DENSITY_CUBE\n")
        if self.bands == True:
            fout.write("\t\t\t&BAND_STRUCTURE\n")
            fout.write("\t\t\t\tADDED_MOS 10\n") 
            fout.write("\t\t\t\tFILE_NAME bands.bs\n")
            fout.write("\t\t\t\t&KPOINT_SET\n")
            fout.write("\t\t\t\t\tUNITS B_VECTOR\n")
            fout.write("\t\t\t\t\tSPECIAL_POINT GAMA 0.0000000000     0.0000000000     0.0000000000\n") # GAMA
            fout.write("\t\t\t\t\tSPECIAL_POINT X    0.5000000000     0.0000000000     0.5000000000\n") # X
            fout.write("\t\t\t\t\tSPECIAL_POINT W    0.5000000000     0.2500000000     0.7500000000\n") # W
            fout.write("\t\t\t\t\tSPECIAL_POINT GAMA 0.0000000000     0.0000000000     0.0000000000\n") # GAMA
            fout.write("\t\t\t\t\tSPECIAL_POINT U    0.6250000000     0.2500000000     0.6250000000\n") # U
            fout.write("\t\t\t\t\tSPECIAL_POINT W    0.5000000000     0.2500000000     0.7500000000\n") # W
            fout.write("\t\t\t\t\tNPOINTS 20\n")
            fout.write("\t\t\t\t&END KPOINT_SET\n")
            fout.write("\t\t\t&END BAND_STRUCTURE\n")
        if self.pdos == True:
            fout.write("\t\t\t&PDOS\n")
            fout.write("\t\t\t\tNLUMO -1\n") # print all projected DOS available
            fout.write("\t\t\t\tCOMPONENTS\n") # split the density by quantum number
            fout.write("\t\t\t&END PDOS\n")
        if self.moments == True:
            fout.write("\t\t\t&MOMENTS\n")
            fout.write("\t\t\t\tPERIODIC %s\n" % False)
            fout.write("\t\t\t&END MOMENTS\n")
        fout.write("\t\t&END PRINT\n")
    def print_electron_density(self):
        self.electron_density = True

    def print_bands(self):
        self.bands = True
    def print_pdos(self):
        self.pdos = True
    def print_moments(self):
        self.moments = True

    def set_params(self, params):
        for item in params:
            if len(item.spit("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass
