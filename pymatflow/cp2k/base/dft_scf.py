#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

# =============================
# CP2K / FORCE_EVAL / DFT / SCF
# =============================

class cp2k_dft_scf_diagonalization_davidson:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
       """
       fout: a file stream for writing
       """
       fout.write("\t\t\t\t&DAVIDSON\n")
       for item in self.params:
           if self.params[item] is not None:
               fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
       fout.write("\t\t\t\t&END DAVIDSON\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_dft_scf_diagonalization_diag_sub_scf_mixing:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
       """
       fout: a file stream for writing
       """
       fout.write("\t\t\t\t\t&MIXING\n")
       for item in self.params:
           if self.params[item] is not None:
               fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
       fout.write("\t\t\t\t\t&END MIXING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_scf_diagonalization_diag_sub_scf:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
       """
       fout: a file stream for writing
       """
       fout.write("\t\t\t\t&DIAG_SUB_SCF\n")
       for item in self.params:
           if self.params[item] is not None:
               fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
       fout.write("\t\t\t\t&END DIAG_SUB_SCF\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "MIXING":
                self.mixing.set_params({item: params[item]})


class cp2k_dft_scf_diagonalization_filter_matrix:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
       """
       fout: a file stream for writing
       """
       fout.write("\t\t\t\t&FILTER_MATRIX\n")
       for item in self.params:
           if self.params[item] is not None:
               fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
       fout.write("\t\t\t\t&END FILTER_MATRIX\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_dft_scf_diagonalization_krylov:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
       """
       fout: a file stream for writing
       """
       fout.write("\t\t\t\t&KRYLOV\n")
       for item in self.params:
           if self.params[item] is not None:
               fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
       fout.write("\t\t\t\t&END KRYLOV\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_dft_scf_diagonalization_ot:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
       """
       fout: a file stream for writing
       """
       fout.write("\t\t\t\t&OT\n")
       for item in self.params:
           if self.params[item] is not None:
               fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
       fout.write("\t\t\t\t&END OT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_dft_scf_diagonalization:
    def __init__(self):
        self.params = {
                "ALGORITHM": "STANDARD",
                "EPS_ADAPT": None,
                "EPS_ITER": None,
                "EPS_JACOBI": None,
                "JACOBI_THRESHOLD": None,
                "MAX_ITER": None,
                }
        self.status = False

        self.section = "TRUE" # TRUE, FALSE, true, false, True, False, actually regardness of uppercase or lowercase

        self.davidson = cp2k_dft_scf_diagonalization_davidson()
        self.diag_sub_scf = cp2k_dft_scf_diagonalization_diag_sub_scf()
        self.filter_matrix = cp2k_dft_scf_diagonalization_filter_matrix()
        self.krylov = cp2k_dft_scf_diagonalization_krylov()
        self.ot = cp2k_dft_scf_diagonalization_ot()

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&DIAGONALIZATION %s\n" % str(self.section))
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.davidson.status == True:
            self.davidson.to_input(fout)
        if self.diag_sub_scf.status == True:
            self.diag_sub_scf.to_input(fout)
        fout.write("\t\t\t&END DIAGONALIZATION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                if item.split("-")[3] == "DIAG_SUB_SCF":
                    self.diag_sub_scf.section = params[item]
                else:
                    self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "DAVIDSON":
                self.davidson.set_params({item: params[item]})
            elif item.split("-")[3] == "DIAG_SUB_SCF":
                self.diag_sub_scf.set_params({item: params[item]})
            elif item.split("-")[3] == "FILTER_MATRIX":
                self.filter_matrix.set_params({item: params[item]})
            elif item.split("-")[3] == "KRYLOV":
                self.krylov.set_params({item: params[item]})
            elif item.split("-")[3] == "OT":
                self.ot.set_params({item: params[item]})



class cp2k_dft_scf_ot:
    def __init__(self):
        self.params = {
                "MINIMIZER": "DIIS",
                "PRECONDITIONER": "FULL_ALL",
                "ENERGY_GAP": 0.001,
                }
        self.status = False

        self.section = "FALSE" # TRUE, FALSE, true, false, True, False, actually regardness of uppercase or lowercase

    def to_input(self, fout):
        fout.write("\t\t\t&OT %s\n" % (str(self.section)))
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END OT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_scf_mixing:
    def __init__(self):
        self.params = {
                "ALPHA": 0.4,
                "BETA": None,
                "METHOD": "BROYDEN_MIXING",
                }
        self.status = False

        self.section = 'TRUE' # TRUE, FALSE, true, false, True, False, actually regardness of uppercase or lowercase

    def to_input(self, fout):
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

class cp2k_dft_scf_mom:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.section = 'TRUE' # TRUE, FALSE, true, false, True, False, actually regardness of uppercase or lowercase

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&MOM %s\n" % str(self.section))
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END MOM\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_scf_outer_scf_cdft_opt:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&CDFT_OPT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END CDFT_OPT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_scf_outer_scf:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.section = 'TRUE' # TRUE, FALSE, true, false, True, False, actually regardness of uppercase or lowercase

        self.cdft_opt = cp2k_dft_scf_outer_scf_cdft_opt()

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&OUTER_SCF %s\n" % str(self.section))
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.cdft_opt.status == True:
            self.cdft_opt.to_input(fout)
        fout.write("\t\t\t&END OUTER_SCF\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "CDFT_OPT":
                self.cdft_opt.set_params({item: params[item]})


class cp2k_dft_scf_smear:
    def __init__(self):
        self.params = {
                "METHOD": None,
                "ELECTRONIC_TEMPERATURE": None,
                "WINDOW_SIZE": None
                }
        self.status = False

        self.section = "FALSE" # TRUE, FALSE, true, false, True, False, actually regardness of uppercase or lowercase

        # basic_setting
        self.section = "FALSE" #False
        self.params["METHOD"] = 'FERMI_DIRAC'
        self.params["ELECTRONIC_TEMPERATURE"] = 300.0


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&SMEAR %s\n" % str(self.section))
        if self.section.upper() == "TRUE":
            if self.params["METHOD"] is not None:
                fout.write("\t\t\t\tMETHOD %s\n" % self.params["METHOD"])
            if self.params["METHOD"] != None and self.params["METHOD"].upper() == "ENERGY_WINDOW" and self.params["WINDOW_SIZE"] != None:
                fout.write("\t\t\t\tWINDOW_SIZE %f\n" % self.params["WINDOW_SIZE"])
            elif self.params["METHOD"] != None and self.params["METHOD"].upper() == "FERMI_DIRAC" and self.params["ELECTRONIC_TEMPERATURE"] != None:
                fout.write("\t\t\t\tELECTRONIC_TEMPERATURE %f\n" % self.params["ELECTRONIC_TEMPERATURE"])
        fout.write("\t\t\t&END SMEAR\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_scf_print_davidson_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_print_davidson:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_scf_print_davidson_each()
        # basic setting
        self.each.status = True
        
        self.section = "HIGH" # DEBUG | HIGH | LOW | MEDIUM | OFF | ON | SILENT
    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&DAVIDSON %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True: 
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END DAVIDSON\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})

class cp2k_dft_scf_print_detailed_energy_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_print_detailed_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_scf_print_detailed_energy_each()
        # basic setting
        self.each.status = True
    
        self.section = "HIGH" # DEBUG | HIGH | LOW | MEDIUM | OFF | ON | SILENT

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&DETAILED_ENERGY %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True: 
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END DETAILED_ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})

class cp2k_dft_scf_print_diag_sub_scf_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_print_diag_sub_scf:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_scf_print_diag_sub_scf_each()
        # basic setting
        self.each.status = True
        
        self.section = "HIGH" # DEBUG | HIGH | LOW | MEDIUM | OFF | ON | SILENT
    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&DIAG_SUB_SCF %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True: 
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END DIAG_SUB_SCF\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})

class cp2k_dft_scf_print_diis_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_print_diis_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_scf_print_diis_info_each()
        # basic setting
        self.each.status = True
        
        self.section = "HIGH" # DEBUG | HIGH | LOW | MEDIUM | OFF | ON | SILENT
    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&DIIS_INFO %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True: 
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END DIIS_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})

class cp2k_dft_scf_print_filter_matrix_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_print_filter_matrix:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_scf_print_filter_matrix_each()
        # basic setting
        self.each.status = True
        
        self.section = "HIGH" # DEBUG | HIGH | LOW | MEDIUM | OFF | ON | SILENT
    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&FILTER_MATRIX %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True: 
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END FILTER_MATRIX\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})



class cp2k_dft_scf_print_iteration_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_print_iteration_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_scf_print_iteration_info_each()
        # basic setting
        self.each.status = True
        
        self.section = "HIGH" # DEBUG | HIGH | LOW | MEDIUM | OFF | ON | SILENT
    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&ITERATION_INFO %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True: 
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END ITERATION_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})


class cp2k_dft_scf_print_lanczos_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_print_lanczos:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_scf_print_lanczos_each()
        # basic setting
        self.each.status = True
        
        self.section = "HIGH" # DEBUG | HIGH | LOW | MEDIUM | OFF | ON | SILENT
    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&LANCZOS %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True: 
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END LANCZOS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})

class cp2k_dft_scf_print_mos_molden_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_print_mos_molden:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_scf_print_mos_molden_each()
        # basic setting
        self.each.status = True
        
        self.section = "HIGH" # DEBUG | HIGH | LOW | MEDIUM | OFF | ON | SILENT
    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&MOS_MOLDEN %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True: 
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END MOS_MOLDEN\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})


class cp2k_dft_scf_print_mo_magnitude_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_print_mo_magnitude:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_scf_print_mo_magnitude_each()
        # basic setting
        self.each.status = True
        
        self.section = "HIGH" # DEBUG | HIGH | LOW | MEDIUM | OFF | ON | SILENT
    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&MO_MAGNITUDE %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True: 
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END MO_MAGNITUDE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})

class cp2k_dft_scf_print_mo_orthonormality_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_print_mo_orthonormality:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_scf_print_mo_orthonormality_each()
        # basic setting
        self.each.status = True
    
        self.section = "HIGHG" # DEBUG | HIGH | LOW | MEDIUM | OFF | ON | SILENT


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&MO_ORTHONORMALITY %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True: 
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END MO_ORTHONORMALITY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})

class cp2k_dft_scf_print_program_run_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_print_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_scf_print_program_run_info_each()
        # basic setting
        self.each.status = True
    
        self.section = "HIGH" # DEBUG | HIGH | LOW | MEDIUM | OFF | ON | SILENT


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PROGRAM_RUN_INFO %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True: 
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})


class cp2k_dft_scf_print_restart_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_print_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_scf_print_restart_each()
        # basic setting
        self.each.status = True
    
        self.section = "HIGH" # DEBUG | HIGH | LOW | MEDIUM | OFF | ON | SILENT


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&RESTART %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True: 
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END RESTART\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})

class cp2k_dft_scf_print_restart_history_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_print_restart_history:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_scf_print_restart_history_each()
        # basic setting
        self.each.status = True
    
        self.section = "HIGH" # DEBUG | HIGH | LOW | MEDIUM | OFF | ON | SILENT

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&RESTART_HISTORY %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True: 
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END RESTART_HISTORY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})


class cp2k_dft_scf_print_total_densities_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_dft_scf_print_total_densities:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_scf_print_total_densities_each()
        # basic setting
        self.each.status = True

        self.section = "HIGH" # DEBUG | HIGH | LOW | MEDIUM | OFF | ON | SILENT
    
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&TOTAL_DENSITIES %s\n" % self.section)
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True: 
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END TOTAL_DENSITIES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})


class cp2k_dft_scf_print:
    def __init__(self):
        self.params = {
                "DM_RESTART_WRITE": None,
                }
        self.status = False

        self.davidson = cp2k_dft_scf_print_davidson()
        self.detailed_energy = cp2k_dft_scf_print_detailed_energy()
        self.diag_sub_scf = cp2k_dft_scf_print_diag_sub_scf()
        self.diis_info = cp2k_dft_scf_print_diis_info()
        self.filter_matrix = cp2k_dft_scf_print_filter_matrix()
        self.iteration_info = cp2k_dft_scf_print_iteration_info()
        self.lanczos = cp2k_dft_scf_print_lanczos()
        self.mos_molden = cp2k_dft_scf_print_mos_molden()
        self.mo_magnitude = cp2k_dft_scf_print_mo_magnitude()
        self.mo_orthonormality = cp2k_dft_scf_print_mo_orthonormality()
        self.program_run_info = cp2k_dft_scf_print_program_run_info()
        self.restart = cp2k_dft_scf_print_restart()
        self.restart_history = cp2k_dft_scf_print_restart_history()
        self.total_densities = cp2k_dft_scf_print_total_densities()

        # basic setting
        self.restart.status = True
        self.restart_history.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.davidson.status == True:
            self.davidson.to_input(fout)
        if self.detailed_energy.status == True:
            self.detailed_energy.to_input(fout)
        if self.diag_sub_scf.status == True:
            self.diag_sub_scf.to_input(fout)
        if self.diis_info.status == True:
            self.diis_info.to_input(fout)
        if self.filter_matrix.status == True:
            self.filter_matrix.to_input(fout)
        if self.iteration_info.statu == True:
            self.iteration_info.to_input(fout)
        if self.lanczos.status == True:
            self.lanczos.to_input(fout)
        if self.mos_molden.status == True:
            self.mos_molden.to_input(fout)
        if self.mo_magnitude.status == True:
            self.mo_magnitude.to_input(fout)
        if self.orthonormality.status == True:
            self.orthonormality.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.restart.status == True:
            self.restart.to_input(fout)
        if self.restart_history.status == True:
            self.restart_history.to_input(fout)
        if self.total_densities.status == True:
            self.total_densities.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params: 
            if len(item.split("-")) == 4:
                if item.split("-")[3] == "DAVIDSON":
                    self.davidson.section = params[item]
                elif item.split("-")[3] == "DETAILED_ENERGY":
                    self.detailed_energy.section = params[item]
                elif item.split("-")[3] == "DIAG_SUB_SCF":
                    self.diag_sub_scf.section = params[item]
                elif item.split("-")[3] == "DIIS_INFO":
                    self.diis_info.section = params[item]
                elif item.split("-")[3] == "FILTER_MATRIX":
                    self.filter_matrix.section = params[item]
                elif item.split("-")[3] == "ITERATION_INFO":
                    sel.fiteration_info.section = params[item]
                elif item.split("-")[3] == "LANCZOS":
                    self.lanczos.section = params[item]
                elif item.split("-")[3] == "MOS_MOLDEN":
                    self.mos_molden.section = params[item]
                elif item.split("-")[3] == "MO_MAGNITUDE":
                    self.mo_magnitude.section = params[item]
                elif item.split("-")[3] == "MO_ORTHONORMALITY":
                    self.mo_orthonormality.section = params[item]
                elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                    self.program_run_info.section = params[item]
                elif item.split("-")[3] == "RESTART":
                    self.restart.section = params[item]
                elif item.split("-")[3] == "RESTART_HISTORY":
                    self.restart_history.section = params[item]
                elif item.split("-")[3] == "TOTAL_DENSITIES":
                    self.total_densities.section = params[item]
                else:
                    self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "DAVIDSON":
                self.dividson.set_params({item: params[item]})
            elif item.split("-")[3] == "DETAILED_ENERGY":
                self.detailed_energy.set_params({item: params[item]})
            elif item.split("-")[3] == "DIAG_SUB_SCF":
                self.diag_sub_scf.set_params({item: params[item]})
            elif item.split("-")[3] == "DIIS_INFO":
                self.diis_info.set_params({item: params[item]})
            elif item.split("-")[3] == "FILTER_MATRIX":
                self.filter_matrix.set_params({item: params[item]})
            elif item.split("-")[3] == "ITERATION_INFO":
                self.iteration_info.set_params({item: params[item]})
            elif item.split("-")[3] == "LANCZOS":
                self.lanczos.set_params({item: params[item]})
            elif item.split("-")[3] == "MOS_MOLDEN":
                self.mos_molden.set_params({item: params[item]})
            elif item.split("-")[3] == "MO_MAGNITUDE":
                self.mo_magnitude.set_params({item: params[item]})
            elif item.split("-")[3] == "MO_ORTHONORMALITY":
                self.mo_orthonormality.set_params({item: params[item]})
            elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[3] == "RESTART":
                self.restart.set_params({item: params[item]})
            elif item.split("-")[3] == "RESTART_HISTORY":
                self.restart_history.set_params({item: params[item]})
            elif item.split("-")[3] == "TOTAL_DENSITIES":
                self.total_densities.set_params({item: params[item]})


class cp2k_dft_scf:
    def __init__(self):
        self.params = {
                "ADDED_MOS": 20,
                "SCF_GUESS": "ATOMIC",
                "EPS_SCF": 1.0E-07,
                "MAX_SCF": 50,
                "MAX_DIIS": None,
                "MAX_ITER_LUMO": None,
                "MAX_SCF_HISTORY": None,
                "ROKS_SCHEME": None,
                }
        self.status = False

        self.diagonalization = cp2k_dft_scf_diagonalization()
        self.ot = cp2k_dft_scf_ot()
        self.mom = cp2k_dft_scf_mom()
        self.mixing = cp2k_dft_scf_mixing()
        self.smear = cp2k_dft_scf_smear()
        self.printout = cp2k_dft_scf_print()
        self.outer_scf = cp2k_dft_scf_outer_scf()

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&SCF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.diagonalization.section.upper() == "TRUE" and self.ot.section.upper() == "FALSE":
            self.diagonalization.to_input(fout)
            self.mixing.to_input(fout)
        elif self.ot.section.upper() == "TRUE" and self.diagonalization.section.upper() == "FALSE":
            self.ot.to_input(fout)
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
        if self.smear.section != None and self.smear.section.upper() == "TRUE":
            """
            if self.params["ADDED_MOS"] == None or self.params["ADDED_MOS"] == 0:
                print("==========================================================\n")
                print("                      Warning !!!\n")
                print("==========================================================\n")
                print("If you are using smearing, you should set ADDED_MOS too!!!\n")
                sys.exit(1)
            """
            self.smear.to_input(fout)
        
        if self.outer_scf.section != None and self.outer_scf.section.upper() == "TRUE":
            self.outer_scf.to_input(fout)

        if self.printout.status == True:
            self.printout.to_input(fout)
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
                elif item.split("-")[-1] == "OUTER_SCF":
                    self.outer_scf.section = params[item]
                elif item.split("-")[-1] == "MOM":
                    self.mom.section = params[item]
                elif item.split("-")[-1] == "MIXING":
                    self.mixing.section = params[item]
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

