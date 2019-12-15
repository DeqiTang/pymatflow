#!/usr/bin/env python
# _*_ coding: utf-8 _*_


# =================================
# CP2K / FORCE_EVAL /DFT / LOCALIZE
# =================================

class cp2k_dft_localize_print_loc_restart_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_localize_print_loc_restart:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_localize_print_loc_restart_each()

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t&LOC_RESTART\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END LOC_RESTART\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_localize_print_molecular_dipoles_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_localize_print_molecular_dipoles:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_localize_print_molecular_dipoles_each()

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t&MOLECULAR_DIPOLES\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END MOLECULAR_DIPOLES\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_localize_print_molecular_states_cubes_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_localize_print_molecular_states_cubes:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_localize_print_molecular_states_cubes_each()

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t\t&CUBES\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END CUBES\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_localize_print_molecular_states_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_localize_print_molecular_states:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cubes = cp2k_dft_localize_print_molecular_states_cubes()
        self.each = cp2k_dft_localize_print_molecular_states_each()

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t&MOLECULAR_STATES\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.cubes.status == True:
            self.cubes.to_input(fout)
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END MOLECULAR_STATES\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CUBES":
                self.cubes.set_params({item: params[item]})
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_localize_print_program_run_info_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_localize_print_program_run_info:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_localize_print_program_run_info_each()

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END PROGRAM_RUN_INFO\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_localize_print_total_dipole_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_localize_print_total_dipole:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_localize_print_total_dipole_each()

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t&TOTAL_DIPOLE\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END TOTAL_DIPOLE\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_localize_print_wannier_centers_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_localize_print_wannier_centers:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_localize_print_wannier_centers_each()

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t&WANNIER_CENTERS\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END WANNIER_CENTERS\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_localize_print_wannier_cubes_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_localize_print_wannier_cubes:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_localize_print_wannier_cubes_each()

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t&WANNIER_CUBES\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END WANNIER_CUBES\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_localize_print_wannier_spreads_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_localize_print_wannier_spreads:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_localize_print_wannier_spreads_each()

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t&WANNIER_SPREADS\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END WANNIER_SPREADS\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_localize_print_wannier_states_cubes_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_localize_print_wannier_states_cubes:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_localize_print_wannier_states_cubes_each()

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t\t&CUBES\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END CUBES\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_localize_print_wannier_states_each:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_localize_print_wannier_states:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cubes = cp2k_dft_localize_print_wannier_states_cubes()
        self.each = cp2k_dft_localize_print_wannier_states_each()

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t&WANNIER_STATES\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.cubes.status == True:
            self.cubes.to_input(fout)
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END WANNIERS_STATES\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CUBES":
                self.cubes.set_params({item: params[item]})
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass



class cp2k_dft_localize_print:
    """
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&print\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&end print\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_localize:
    """
    about:
        to calculate ir spectra from md running
        it is necessary to have dipole information
        for the molecules available in the simulated
        trajectory.
    """
    def __init__(self):
        self.params = {
                "method": None,
                "max_iter": None,
                }
        self.status = False

        self.printout = cp2k_dft_localize_print()

        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&localize %s\n" % ".true.")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t%s %s\n" % (item, str(self.params[item])))
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t&END LOCALIZE\n")
    
    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass
