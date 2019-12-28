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


class cp2k_optimize_basis_fit_kind_constrain_exponents:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&CONSTRAIN_EXPONENTS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END CONSTRAIN_EXPONENTS\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_optimize_basis_fit_kind_derived_basis_sets:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&DERIVED_BASIS_SETS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END DERIVED_BASIS_SETS\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_optimize_basis_fit_kind:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.constrain_exponents = cp2k_optimize_basis_fit_kind_constrain_exponents()
        self.derived_basis_sets = cp2k_optimize_fit_kind_derived_basis_sets()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&FIT_KIND\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.constrain_exponents.status == True:
            self.constrain_exponents.to_input(fout)
        if self.derived_basis_sets.status == True:
            self.derived_basis_sets.to_input(fout)
        fout.write("\t&END FIT_KIND\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "CONSTRAIN_EXPONENTS":
                self.constrain_exponents.set_params({item: params[item]})
            elif item.split("-")[1] == "DERIVED_BASIS_SETS":
                self.derived_basis_sets.set_params({item: params[item]})
            else:
                pass


class cp2k_optimize_basis_optimization:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&OPTIMIZATION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END OPTIMIZATION\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_optimize_basis_training_files:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&TRAINING_FILES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END TRAINING_FILES\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_optimize_basis:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False
            
        self.fit_kind = cp2k_optimize_basis_fit_kind()
        self.optimization = cp2k_optimize_basis_optimization()
        self.traning_files = cp2k_optimize_basis_training_files()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("&OPTIMIZE_BASIS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t%s %s\n" % (item, self.params[item]))
        if self.fit_kind.status == True:
            self.fit_kind.to_input(fout)
        if self.optimization.status == True:
            self.optimization.to_input(fout)
        if self.training_files.status == True:
            self.training_files.to_input(fout)
        fout.write("&END OPTIMIZE_BASIS\n")
        fout.write("\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 1:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[0] == "FIT_KIND":
                self.fit_kind.set_params({item: params[item]})
            elif item.split("-")[0] == "OPTIMIZATION":
                self.optimization.set_params({item: params[item]})
            elif item.split("-")[0] == "TRAINING_FILES":
                self.training_files.set_params({item: params[item]})
            else:
                pass


