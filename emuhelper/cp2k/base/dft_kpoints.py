#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# =================================
# CP2K / FORCE_EVAL / DFT / KPOINTS
# =================================
class cp2k_dft_kpoints:
    """
    """
    def __init__(self):
        self.params = {
                "EPS_GEO": None,
                "FULL_GRID": None,
                "KPOINT": None,
                "PARALLEL_GROUP_SIZE": None,
                "SCHEME": None,
                "SYMMETRY": None,
                "UNITS": None,
                "VERBOSE": None,
                "WAVEFUNCTIONS": None,
                }
        self.kpoints_mp = [1, 1, 1]
        self.basic_setting()

    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&KPOINTS\n")
        if self.params["SCHEME"] == "MONKHORST-PACK":
            fout.write("\t\t\tSCHEME MONKHORST-PACK %d %d %d\n" % (self.kpoints_mp[0], self.kpoints_mp[1], self.kpoints_mp[2]))
        elif self.params["SCHEME"] == "GAMMA":
            fout.write("\t\t\tSCHEME GAMMA\n")
        fout.write("\t\t&END KPOINTS\n")

    def set_kpoints(self, scheme="MONKHORST-PACK", kpoints_mp=[1, 1, 1]):
        self.kpoints_mp = kpoints_mp
        self.params["SCHEME"] = scheme

    def basic_setting(self):
        self.params["SCHEME"] = "GAMMA"

