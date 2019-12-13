#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys


# =================================
# CP2K / FORCE_EVAL / DFT / KPOINTS
# =================================
class cp2k_dft_kpoints:
    """
    Note:
        the setting of SCHEME should obey the regulation of cp2k.
        like you can set it to be GAMMA, or NONE.
        but when you set it to be MONKHORST-PACK, you should also
        add the three integers behind it, like SCHEME is a string
        as: 'MONKHORST-PACK 3 3 3'.
        we add a check to guaranee this

        OT not possible with kpoint calculations, so when doing
        OT calculation, we must use DFT-KPOINTS-SCHEME = NONE
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
        self.status = False

        # basic_setting
        self.params["SCHEME"] = "GAMMA"

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        # check the kpoints setting, especially the SCHEME
        self.check_kpoints()
        #
        fout.write("\t\t&KPOINTS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END KPOINTS\n")

    def check_kpoints(self):
        if self.params["SCHEME"][0:4].upper() == "MONK":
            if len(self.params["SCHEME"].split()) != 4:
                print("=============================================================\n")
                print("                        WARNING !!!\n")
                print("==============================================================\n")
                print("dft.kpoints.check_kpoints():\n")
                print("when you setting the KPOINTS-SHCEME to\n")
                print("MONKHORST-PACK type, the number of k points\n")
                print("in all 3 dimensions has to be supplied along with the keyword.\n")
                print("E.g. MONKHORST-PACK 3 3 3.")
                sys.exit(1)


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass
