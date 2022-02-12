#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

"""
Usage:
"""

class cp2k_vibrational_analysis_mode_selective_involved_atoms:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&INVOLVED_ATOMS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END INVOLVED_ATOMS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_vibrational_analysis_mode_selective_print_ms_restart_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_vibrational_analysis_mode_selective_print_ms_restart:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_vibrational_analysis_mode_selective_print_ms_restart_each()

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&MS_RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END MS_RESTART\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_vibrational_analysis_mode_selective_print:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.ms_restart = cp2k_vibrational_analysis_mode_selective_print_ms_restart()

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.ms_restart.status == True:
            self.ms_restart.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "MS_RESTART":
                self.ms_restart.set_params({item: params[item]})
            else:
                pass


class cp2k_vibrational_analysis_mode_selective:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.involved_atoms = cp2k_vibrational_analysis_mode_selective_involved_atoms()
        self.printout = cp2k_vibrational_analysis_mode_selective_print()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&MODE_SELECTIVE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.involved_atoms.status == True:
            self.involved_atoms.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t&END MODE_SELECTIVE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "INVOLVED_ATOMS":
                self.involved_atoms.set_params({item: params[item]})
            elif item.split("-")[1] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass


class cp2k_vibrational_analysis_print_banner_each:
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
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_vibrational_analysis_print_banner:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_vibrational_analysis_print_banner_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&BANNER\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END BANNER\n")
        fout.write("\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_vibrational_analysis_print_cartesian_eigs_each:
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
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_vibrational_analysis_print_cartesian_eigs:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_vibrational_analysis_print_cartesian_eigs_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&CARTESIAN_EIGS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END CARTESIAN_EIGS\n")
        fout.write("\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_vibrational_analysis_print_hessian_each:
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
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_vibrational_analysis_print_hessian:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_vibrational_analysis_print_hessian_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&HESSIAN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END HESSIAN\n")
        fout.write("\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_vibrational_analysis_print_molden_vib_each:
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
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_vibrational_analysis_print_molden_vib:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_vibrational_analysis_print_molden_vib_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&MOLDEN_VIB\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END MOLDEN_VIB\n")
        fout.write("\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_vibrational_analysis_print_program_run_info_each:
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
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_vibrational_analysis_print_program_run_info:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_vibrational_analysis_print_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PROGRAM_RUN_INFO\n")
        fout.write("\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_vibrational_analysis_print_rotational_info_each:
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
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_vibrational_analysis_print_rotational_info:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_vibrational_analysis_print_rotational_info_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&ROTATIONAL_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END ROTATIONAL_INFO\n")
        fout.write("\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_vibrational_analysis_print:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.banner = cp2k_vibrational_analysis_print_banner()
        self.cartesian_eigs = cp2k_vibrational_analysis_print_cartesian_eigs()
        self.hessian = cp2k_vibrational_analysis_print_hessian()
        self.molden_vib = cp2k_vibrational_analysis_print_molden_vib()
        self.program_run_info = cp2k_vibrational_analysis_print_program_run_info()
        self.rotational_info = cp2k_vibrational_analysis_print_rotational_info()

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.banner.status == True:
            self.banner.to_input(fout)
        if self.cartesian_eigs.status == True:
            self.cartesian_eigs.to_input(fout)
        if self.hessian.status == True:
            self.hessian.to_input(fout)
        if self.molden_vib.status == True:
            self.molden_vib.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.rotational_info.status == True:
            self.rotational_info.to_input(fout)
        fout.write("\t&END PRINT\n")
        fout.write("\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "BANNER":
                self.banner.set_params({item: params[item]})
            elif item.split("-")[1] == "CARTESIAN_EIGS":
                self.cartesian_eigs.set_params({item: params[item]})
            elif item.split("-")[1] == "HESSIAN":
                self.hessian.set_params({item: params[item]})
            elif item.split("-")[1] == "MOLDEN_VIB":
                self.molden_vib.set_params({item: params[item]})
            elif item.split("-")[1] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[1] == "ROTATIONAL_INFO":
                self.rotational_info.set_params({item: params[item]})
            else:
                pass

class cp2k_vibrational_analysis:
    """
    Note:
        if you are calculating IR, you should print out
        dipole moment throught DFT/PRINT/MOMENTS.
        for direct diagonalization method dipole moment
        calculation of Berry phase approach is not supported
        , and for OT method, we can use PERIODIC .TRUE in
        DFT/PRINT/MOMENTS

        inf DFT/PRINT/MOMENTS
        PERIODIC:
            Use Berry phase formula (PERIODIC=T) or simple operator (PERIODIC=F). 
            The latter normally requires that the CELL is periodic NONE.
        so if we want to use Berry phase to calculate dipole moment 
        we have to use OT in SCF, and if we use PERIODIC=F, we have to set 
        PERIODIC=NONE in CELL.
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.mode_selective = cp2k_vibrational_analysis_mode_selective()
        self.printout = cp2k_vibrational_analysis_print()
        # basic setting
        self.params["INTENSITIES"] = True
        self.params["DX"] = 1.0e-2 # default value in Bohr
        self.params["TC_PRESSURE"] = 1.01325e5  # Pa
        self.params["TC_TEMPERATURE"] = 2.7315e2 # K
        self.params["THERMOCHEMISTRY"] = False # Calculation of the thermochemical data. Valid for molecules in the gas phase.


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("&VIBRATIONAL_ANALYSIS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t%s %s\n" % (item, self.params[item]))
        if self.mode_selective.status == True:
            self.mode_selective.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("&END VIBRATIONAL_ANALYSIS\n")
        fout.write("\n")


    def set_params(self, params):
        """
        parameters for sub section(like), are handled over 
        to sub section controllers.
        """
        for item in params:
            if len(item.split("-")) == 1:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[0] == "MODE_SELECTIVE":
                self.mode_selection.set_params({item: params[item]})
            elif item.split("-")[0] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass
