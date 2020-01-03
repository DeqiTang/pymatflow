#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


"""
Usage:
"""



class cp2k_test_cp_dbcsr:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&CP_DBCSR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END CP_DBCSR\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_test_cp_fm_gemm:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&CP_FM_GEMM\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END CP_FM_GEMM\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_test_eigensolver:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&EIGENSOLVER\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END EIGENSOLVER\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_test_eri_mme_test_eri_mme_cutoff_calib:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&CUTOFF_CALIB\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END CUTOFF_CALIB\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_test_eri_mme_test_eri_mme_eri_mme_info_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_test_eri_mme_test_eri_mme_eri_mme_info:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_test_eri_mme_test_eri_mme_eri_mme_info_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&ERI_MME_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END ERI_MME_INFO\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_test_eri_mme_test_eri_mme:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cutoff_calib = cp2k_test_eri_mme_test_eri_mme_cutoff_calib()
        self.eri_mme_info = cp2k_test_eri_mme_test_eri_mme_eri_mme_info()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&ERI_MME\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.cutoff_calib.status == True:
            self.cutoff_calib.to_input(fout)
        if self.eri_mme_info.status == True:
            self.eri_mme_info.to_input(fout)
        fout.write("\t\t&END ERI_MME\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "CUTOFF_CALIB":
                self.cutoff_calib.set_params({item: params[item]})
            elif item.split("-")[2] == "ERI_MME_INFO":
                self.eri_mme_info.set_params({item: params[item]})
            else:
                pass


class cp2k_test_eri_mme_test:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.eri_mme = cp2k_test_eri_mme_test_eri_mme()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&ERI_MME_TEST\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.eri_mme.status == True:
            self.eri_mme.to_input(fout)
        fout.write("\t&END ERI_MME_TESTT\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "ERI_MME":
                self.eri_mme.set_params({item: params[item]})
            else:
                pass


class cp2k_test_grid_information_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_test_grid_information:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.each = cp2k_test_grid_information_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&GRID_INFORMATION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t&END GRID_INFORMATION\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_test_program_run_info_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_test_program_run_info:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_test_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fotu)
        fout.write("\t&END PROGRAM_RUN_INFO\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_test_pw_transfer:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&PW_TRANSFER\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END PW_TRANSFER\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_test_rs_pw_transfer_rs_grid:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&RS_GRID\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END RS_GRID\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_test_rs_pw_transfer:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.rs_grid = cp2k_test_rs_pw_transfer_rs_grid()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&RS_PW_TRANSFER\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.rs_grid.status == True:
            self.rs_grid.to_input(fout)
        fout.write("\t&END RS_PW_TRANSFER\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "RS_GRID":
                self.rs_grid.set_params({item: params[item]})
            else:
                pass


class cp2k_test_shg_integrals_test_basis:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&BASIS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END BASIS\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_test_shg_integrals_test:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.basis = cp2k_test_shg_integrals_test_basis()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&SHG_INTEGRALS_TEST\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.basis.status == True:
            self.basis.to_input(fout)
        fout.write("\t&END SHG_INTEGRALS_TEST\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "BASIS":
                self.basis.set_params({item: params[item]})
            else:
                pass



class cp2k_test:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cp_dbcsr = cp2k_test_cp_dbcsr()
        self.cp_fm_gemm = cp2k_test_cp_fm_gemm()
        self.eigensolver = cp2k_test_eigensolver()
        self.eri_mme_test = cp2k_test_eri_mme_test()
        self.grid_information = cp2k_test_grid_information()
        self.program_run_info = cp2k_test_program_run_info()
        self.pw_transfer = cp2k_test_pw_transfer()
        self.rs_pw_transfer = cp2k_test_rs_pw_transfer()
        self.shg_integrals_test = cp2k_test_shg_integrals_test()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("&TEST\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t%s %s\n" % (item, self.params[item]))
        if self.cp_dbcsr.status == True:
            slef.cp_dbcsr.to_input(fout)
        if self.cp_fm_gemm.status == True:
            self.cp_fm_gemm.to_input(fout)
        if self.eigensolver.status == True:
            self.eigensolver.to_input(fout)
        if self.eri_mme_test.status == True:
            self.eri_mme_test.to_input(fout)
        if self.grid_information.status == True:
            self.grid_information.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.pw_transfer.status == True:
            self.pw_transfer.to_input(fout)
        if self.rs_pw_transfer.status == True:
            self.rs_pw_transfer.to_input(fout)
        if self.shg_integrals_test.status == True:
            self.shg_integrals_test.to_input(fout)
        fout.write("&END TEST\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 1:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[0] == "CP_DBCSR":
                self.cp_dbcsr.set_params({item: params[item]})
            elif item.split("-")[0] == "CP_FM_GEMM":
                self.cp_fm_gemm.set_params({item: params[item]})
            elif item.split("-")[0] == "EIGENSOVLER":
                self.eigensolver.set_params({item: params[item]})
            elif item.split("-")[0] == "ERI_MME_TEST":
                self.eri_mme_test.set_params({item: parmas[item]})
            elif item.split("-")[0] == "GRID_INFORMATION":
                self.grid_information.set_params({item: parmas[item]})
            elif item.split("-")[0] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[0] == "PW_TRANSFER":
                self.pw_transfer.set_params({item: params[item]})
            elif item.split("-")[0] == "RS_PW_TRANSFER":
                self.rs_pw_transfer.set_params({item: parmas[item]})
            elif item.split("-")[0] == "SHG_INTEGRALS_TEST":
                self.shg_integrals_test.set_params({item: params[item]})
            else:
                pass


