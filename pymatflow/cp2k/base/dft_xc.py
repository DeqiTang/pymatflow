#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# =============================
# CP2K / FORCE_EVAL / DFT / XC
# =============================

class cp2k_dft_xc_vdw_potential_pair_potential:
    def __init__(self):
        self.params = {
                }

        self.params["PARAMETER_FILE_NAME"] = "dftd3.dat"

    def to_vdw(self, fout):
        fout.write("\t\t\t\t&PAIR_POTENTIAL\n")
        for item in self.params:
            fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END PAIR_POTENTIAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_vdw_potential_non_local:
    def __init__(self):
        self.params = {
                }

    def to_vdw(self, fout):
        fout.write("\t\t\t\t&NON_LOCAL\n")
        for item in self.params:
            fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END NON_LOCAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_dft_xc_vdw_potential:
    def __init__(self):
        self.section = "FALSE"
        self.params = {
                "POTENTIAL_TYPE": None,
                }
        self.params["POTENTIAL_TYPE"] = "PAIR_POTENTIAL"
        self.pair_potential = cp2k_dft_xc_vdw_potential_pair_potential()
        self.non_local = cp2k_dft_xc_vdw_potential_non_local()

    def to_xc(self, fout):
        fout.write("\t\t\t&VDW_POTENTIAL\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.params["POTENTIAL_TYPE"].upper() == "PAIR_POTENTIAL":
            self.pair_potential.to_vdw(fout)
        elif self.params["POTENTIAL_TYPE"].upper() == "NON_LOCAL":
            self.non_local.to_vdw(fout)
        fout.write("\t\t\t&END VDW_POTENTIAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "PAIR_POTENTIAL":
                self.pair_potential.set_params({item: params[item]})
            elif item.split("-")[3] == "NON_LOCAL":
                self.non_local.set_params({item: params[item]})

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
        self.vdw_potential = cp2k_dft_xc_vdw_potential()

    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&XC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s" % (item, str(self.params[item])))
        self.xc_functional.to_xc(fout)
        if self.vdw_potential.section.upper() == "TRUE":
            self.vdw_potential.to_xc(fout)
        fout.write("\t\t&END XC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                if item.split("-")[-1] == "XC_FUNCTIONAL":
                    self.xc_functional.section = params[item]
                elif item.split("-")[-1] == "VDW_POTENTIAL":
                    self.vdw_potential.section = params[item]
                else:
                    self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "XC_FUNCTIONAL" and len(item.split("-")) > 3:
                self.xc_functional.set_params({item: params[item]})
            elif item.split("-")[2] == "VDW_POTENTIAL":
                self.vdw_potential.set_params({item: params[item]})

