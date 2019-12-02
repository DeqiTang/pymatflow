#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

# =============================
# CP2K / FORCE_EVAL / DFT / SCF
# =============================
class cp2k_dft_scf_diagonalization:
    def __init__(self):
        self.section = "TRUE" # TRUE, FALSE, true, false, True, False, actually regardness of uppercase or lowercase
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
       fout.write("\t\t\t&DIAGONALIZATION %s\n" % str(self.section))
       for item in self.params:
           if self.params[item] is not None:
               fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
       fout.write("\t\t\t&END DIAGONALIZATION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_ot:
    def __init__(self):
        self.section = "FALSE" # TRUE, FALSE, true, false, True, False, actually regardness of uppercase or lowercase
        self.params = {
                "MINIMIZER": "DIIS",
                "PRECONDITIONER": "FULL_ALL",
                "ENERGY_GAP": 0.001,
                }

    def to_scf(self, fout):
        fout.write("\t\t\t&OT %s\n" % (str(self.section)))
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END OT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_scf_mixing:
    def __init__(self):
        self.section = 'TRUE' # TRUE, FALSE, true, false, True, False, actually regardness of uppercase or lowercase
        self.params = {
                "ALPHA": 0.4,
                "BETA": None,
                "METHOD": "BROYDEN_MIXING",
                }
    def to_scf(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&MIXING %s\n" % str(self.section))
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
        self.section = "FALSE" # TRUE, FALSE, true, false, True, False, actually regardness of uppercase or lowercase
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
        if self.section.upper() == "TURE":
            fout.write("\t\t\t\tMETHOD %s\n" % self.params["METHOD"])
            if self.params["METHOD"] == "ENERGY_WINDOW":
                fout.write("\t\t\t\tWINDOW_SIZE %f\n" % self.params["WINDOW_SIZE"])
            elif self.params["METHOD"] == "FERMI_DIRAC":
                fout.write("\t\t\t\tELECTRONIC_TEMPERATURE %f\n" % self.params["ELECTRONIC_TEMPERATURE"])
        fout.write("\t\t\t&END SMEAR\n")

    def basic_setting(self):
        self.section = "FALSE" #False
        self.params["METHOD"] = 'FERMI_DIRAC'
        self.params["ELECTRONIC_TEMPERATURE"] = 300.0

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
        self.ot = cp2k_dft_scf_ot()
        self.mixing = cp2k_dft_scf_mixing()
        self.smear = cp2k_dft_scf_smear()
        self.printout = cp2k_dft_scf_print()

    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&SCF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.diagonalization.section.upper() == "TRUE" and self.ot.section.upper() == "FALSE":
            self.diagonalization.to_scf(fout)
            self.mixing.to_scf(fout)
        elif self.ot.section.upper() == "TRUE" and self.diagonalization.section.upper() == "FALSE":
            self.ot.to_scf(fout)
            #self.mixing.to_scf(fout)
        else:
            print("======================================\n")
            print("             WARNING !!!\n")
            print("your setting:\n")
            print("DFT/SCF/DIAGONALIZATION %s\n" % str(self.diagonalization.section))
            print("DFT/SCF/OT %s\n" % str(self.ot.section))
            print("you can choose either DIAGONALIZATION\n")
            print("or ORBITAL TRANSFORMATION for SCF\n")
            print("======================================\n")
            sys.exit(1)
        if self.smear.section.upper() == "TRUE":
            if self.params["ADDED_MOS"] == None or self.params["ADDED_MOS"] == 0:
                print("If you are using smearing, you should set ADDED_MOS too!!!\n")
                sys.exit(1)
            self.smear.to_scf(fout)
        self.printout.to_scf(fout)
        fout.write("\t\t&END SCF\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 3:
                if item.split("-")[-1] == "SMEAR":
                    self.smear.section = params[item]
                elif item.split("-")[-1] == "DIAGONALIZATION":
                    self.diagonalization.section = params[item]
                elif item.split("-")[-1] == "OT":
                    self.ot.section = params[item]
                else:
                    self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "DIAGONALIZATION":
                self.diagonalization.set_params({item: params[item]})
            elif item.split("-")[2] == "OT":
                self.ot.set_params({item: params[item]})
            elif item.split("-")[2] == "MIXING":
                self.mixing.set_params({item: params[item]})
            elif item.split("-")[2] == "SMEAR" and len(item.split("-")) > 3:
                self.smear.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})

