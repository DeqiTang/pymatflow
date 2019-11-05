#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

# =============================
# CP2K / FORCE_EVAL / DFT / SCF
# =============================
class cp2k_dft_scf_diagonalization:
    def __init__(self):
        self.section = '.TRUE.'
        self.params = {
                "ALGORITHM": "STANDARD",
                "EPS_ADAPT": None,
                "EPS_ITER": None,
                "EPS_JACOBI": None,
                "JACOBI_THRESHOLD": None,
                "MAX_ITER": None,
                }
    def to_scf(self, fout):
       """
       fout: a file stream for writing
       """
       fout.write("\t\t\t&DIAGONALIZATION\n")
       for item in self.params:
           if self.params[item] is not None:
               fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
       fout.write("\t\t\t&END DIAGONALIZATION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_dft_scf_mixing:
    def __init__(self):
        self.section = '.TRUE.'
        self.params = {
                "ALPHA": 0.4,
                "BETA": None,
                "METHOD": "BROYDEN_MIXING",
                }
    def to_scf(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&MIXING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END MIXING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_dft_scf_smear:
    def __init__(self):
        self.section = True # True or False or '.TRUE.' or '.FALSE.'
        self.params = {
                "METHOD": None,
                "ELECTRONIC_TEMPERATURE": None,
                "WINDOW_SIZE": None
                }
        self.basic_setting()

    def to_scf(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&SMEAR %s\n" % str(self.section))
        if self.section == True or self.section == '.TRUE.':
            fout.write("\t\t\t\tMETHOD %s\n" % self.params["METHOD"])
            if self.params["METHOD"] == "ENERGY_WINDOW":
                fout.write("\t\t\t\tWINDOW_SIZE %f\n" % self.params["WINDOW_SIZE"])
            elif self.params["METHOD"] == "FERMI_DIRAC":
                fout.write("\t\t\t\tELECTRONIC_TEMPERATURE %f\n" % self.params["ELECTRONIC_TEMPERATURE"])
        fout.write("\t\t\t&END SMEAR\n")

    def basic_setting(self):
        self.section = True
        self.params["METHOD"] = 'FERMI_DIRAC'
        self.params["ELECTRONIC_TEMPERATURE"] = 300.0

    def set_smear(self, method="FERMI_DIRAC", electronic_temperature=300.0, window_size=0.0):
        self.params["METHOD"] = method
        self.params["ELECTRONIC_TEMPERATURE"] = electronic_temperature
        self.params["WINDOW_SIZE"] = window_size

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_scf_print:
    def __init__(self):
        self.params = {
                "DM_RESTART_WRITE": None,
                }
    def to_scf(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_scf:
    def __init__(self):
        self.params = {
                "ADDED_MOS": 20,
                "SCF_GUESS": "ATOMIC",
                "EPS_SCF": "1.0E-06",
                "MAX_SCF": 50,
                "MAX_DIIS": None,
                "MAX_ITER_LUMO": None,
                "MAX_SCF_HISTORY": None,
                "ROKS_SCHEME": None,
                }
        self.diagonalization = cp2k_dft_scf_diagonalization()
        self.mixing = cp2k_dft_scf_mixing()
        self.smear = cp2k_dft_scf_smear()
        self.printout = cp2k_dft_scf_print()
        self.ifsmear = True

    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&SCF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        self.diagonalization.to_scf(fout)
        self.mixing.to_scf(fout)
        if self.ifsmear == True:
            if self.params["ADDED_MOS"] == None or self.params["ADDED_MOS"] == 0:
                print("If you are using smearing, you should set ADDED_MOS too!!!\n")
                sys.exit(1)
            self.smear.to_scf(fout)
        self.printout.to_scf(fout)
        fout.write("\t\t&END SCF\n")

    def set_params(self, params):
        #self.diagonalization.set_params(params)
        #self.smear.set_params(params)
        #self.mixing.set_params(params)
        for item in params:
            if len(item.split("-")) == 3:
                if item.split("-")[-1] == "SMEAR":
                    self.smear.section = params[item]
                else:
                    self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "DIAGONALIZATION":
                self.diagonalization.set_params({item: params[item]})
            elif item.split("-")[2] == "MIXING":
                self.mixing.set_params({item: params[item]})
            elif item.split("-")[2] == "SMEAR" and len(item.split("-")) > 3:
                self.smear.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})

