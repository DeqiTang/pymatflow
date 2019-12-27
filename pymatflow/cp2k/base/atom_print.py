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

class cp2k_atom_print_admm_admm_basis_basis:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t&BASIS\n")
        for item in self.params:
            fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END BASIS\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_atom_print_admm_admm_basis:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.basis = cp2k_atom_print_admm_admm_basis_basis()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&ADMM_BASIS\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.basis.status == True:
            self.basis.to_input(fout)
        fout.write("\t\t\t&END ADMM_BASIS\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "BASIS":
                self.basis.set_params({item: params[item]})
            else:
                pass


class cp2k_atom_print_admm_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_admm:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_admm_each()
        self.admm_basis = cp2k_atom_print_admm_admm_basis()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&ADMM\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        if self.admm_basis.status == True:
            self.admm_basis.to_input(fout)
        fout.write("\t\t&END ADMM\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            elif item.split("-")[2] == "ADMM_BASIS":
                self.admm_basis.set_params({item: params[item]})
            else:
                pass

class cp2k_atom_print_analyze_basis_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_analyze_basis:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_analyze_basis_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&ANALYZE_BASIS\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END ANALYZE_BASIS\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_atom_print_basis_set_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_basis_set:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_basis_set_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&BASIS_SET\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END BASIS_SET\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_atom_print_fit_basis_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_fit_basis:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_fit_basis_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&FIT_BASIS\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END FIT_BASIS\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_atom_print_fit_density_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_fit_density:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_fit_density_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&FIT_DENSITY\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END FIT_DENSITY\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_atom_print_fit_kgpot_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_fit_kgpot:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_fit_kgpot_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&FIT_KGPOT\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END FIT_KGPOT\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_atom_print_fit_pseudo_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_fit_pseudo:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_fit_pseudo_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&FIT_PSEUDO\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END FIT_PSEUDO\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_atom_print_geometrical_response_basis_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_geometrical_response_basis:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_geometrical_response_basis_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&GEOMETRICAL_RESPONSE_BASIS\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END GEOMETRICAL_RESPONSE_BASIS\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_atom_print_method_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_method_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_method_info_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&METHOD_INFO\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END METHOD_INFO\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_atom_print_orbitals_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_orbitals:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_orbitals_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&ORBITALS\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END ORBITALS\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_atom_print_potential_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_potential:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_potential_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&POTENTIAL\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END POTENTIAL\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_atom_print_program_banner_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_program_banner:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_program_banner_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&PROGRAM_BANNER\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PROGRAM_BANNER\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_atom_print_response_basis_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_response_basis:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_response_basis_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&RESPONSE_BASIS\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END RESPONSE_BASIS\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_atom_print_scf_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_scf_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_scf_info_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&SCF_INFO\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END SCF_INFO\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_atom_print_separable_gaussian_pseudo_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_separable_gaussian_pseudo:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_separable_gaussian_pseudo_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&SEPARABLE_GAUSSIAN_PSEUDO\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END SEPARABLE_GAUSSIAN_PSEUDO\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_atom_print_upf_file_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_print_upf_file:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_atom_print_upf_file_each()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&UPF_FILE\n")
        for item in self.params:
            fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END UPF_FILE\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass



class cp2k_atom_print:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.admm = cp2k_atom_print_admm()
        self.analyze_basis = cp2k_atom_print_analyze_basis()
        self.basis_set = cp2k_atom_print_basis_set()
        self.fit_basis = cp2k_atom_print_fit_basis()
        self.fit_density = cp2k_atom_print_fit_density()
        self.fit_kgpot = cp2k_atom_print_fit_kgpot()
        self.fit_pseudo = cp2k_atom_print_fit_pseudo()
        self.geometrical_response_basis = cp2k_atom_print_geometrical_response_basis()
        self.method_info = cp2k_atom_print_method_info()
        self.orbitals = cp2k_atom_print_orbitals()
        self.potential = cp2k_atom_print_potential()
        self.program_banner = cp2k_atom_print_program_banner()
        self.response_basis = cp2k_atom_print_response_basis()
        self.scf_info = cp2k_atom_print_scf_info()
        self.separable_gaussian_pseudo = cp2k_atom_print_separable_gaussian_pseudo()
        self.upf_file = cp2k_atom_print_upf_file()
        # basic setting


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&PRINT\n")
        for item in self.params:
            fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.admm.status == True:
            self.admm.to_input(fout)
        if self.analyze_basis.status == True:
            self.analyze_basis.to_input(fout)
        if self.basis_set.status == True:
            self.basis_set.to_input(fout)
        if self.fit_basis.status == True:
            self.fit_basis.to_input(fout)
        if self.fit_density.status == True:
            self.fit_density.to_input(fout)
        if self.fit_kgpot.status == True:
            self.fit_kgpot.to_input(fout)
        if self.fit_pseudo.status == True:
            self.fit_pseudo.to_input(fout)
        if self.geometrical_response_basis.status == True:
            self.geometrical_response_basis.to_input(fout)
        if self.method_info.status == True:
            self.method_info.to_input(fout)
        if self.orbitals.status == True:
            self.orbitals.to_input(fout)
        if self.potential.status == True:
            self.potential.to_input(fout)
        if self.program_banner.status == True:
            self.program_banner.to_input(fout)
        if self.response_basis.status == True:
            self.response_basis.to_input(fout)
        if self.scf_info.status == True:
            self.scf_info.to_input(fout)
        if self.separable_gaussian_pseudo.status == True:
            self.separable_gaussian_pseudo.to_input(fout)
        if self.upf_file.status == True:
            self.upf_file.to_input(fout)
        fout.write("\t&END PRINT\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "ADMM":
                self.admm.set_params({item: params[item]})
            elif item.split("-")[1] == "ANALYZE_BASIS":
                self.analyze_basis.set_params({item: params[item]})
            elif item.split("-")[1] == "BASIS_SET":
                self.basis_set.set_params({item: params[item]})
            elif item.split("-")[1] == "FIT_BASIS":
                self.fit_basis.set_params({item: params[item]})
            elif item.split("-")[1] == "FIT_DENSITY":
                self.fit_density.set_params({item: params[item]})
            elif item.split("-")[1] == "FIT_KGPOT":
                self.kgpot.set_params({item: params[item]})
            elif item.split("-")[1] == "FIT_PSEUDO":
                self.fit_pseudo.set_params({item: params[item]})
            elif item.split("-")[1] == "GEOMETRICAL_RESPONSE_BASIS":
                self.geometrical_response_basis.set_params({item: params[item]})
            elif item.split("-")[1] == "METHOD_INFO":
                self.method_info.set_params({item: params[item]})
            elif item.split("-")[1] == "ORBITALS":
                self.orbitals.set_params({item: params[item]})
            elif item.split("-")[1] == "POTENTIAL":
                self.potential.set_params({item: params[item]})
            elif item.split("-")[1] == "PROGRAM_BANNER":
                self.program_banner.set_params({item: params[item]})
            elif item.split("-")[1] == "RESPONSE_BASIS":
                self.response_basis.set_params({item: params[item]})
            elif item.split("-")[1] == "SCF_INFO":
                self.scf_info.set_params({item: params[item]})
            elif item.split("-")[1] == "SEPARABLE_GAUSSIAN_PSEUDO":
                self.separable_gaussian_pseudo.set_params({item: params[item]})
            elif item.split("-")[1] == "UPF_FILE":
                self.upf_file.set_params({item: params[item]})
            else:
                pass



