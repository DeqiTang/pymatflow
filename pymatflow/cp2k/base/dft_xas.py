#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ==============================
# ==============================

class cp2k_dft_xas_localize_print_loc_restart_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_localize_print_loc_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_localize_print_loc_restart_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&LOC_RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END LOC_RESTART\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xas_localize_print_molecular_dipoles_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_localize_print_molecular_dipoles:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_localize_print_molecular_dipoles_each()

        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&MOLECULAR_DIPOLES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END MOLECULAR_DIPOLES\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_xas_localize_print_molecular_states_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_localize_print_molecular_states_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_localize_print_molecular_states_cubes_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t&END CUBES\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[6] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_xas_localize_print_molecular_states_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_localize_print_molecular_states:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cubes = cp2k_dft_xas_localize_print_molecular_states_cubes()
        self.each = cp2k_dft_xas_localize_print_molecular_states_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&MOLECULAR_STATES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.cubes.status == True:
            self.cubes.to_input(fout)
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END MOLECULAR_STATES\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "CUBES":
                self.cubes.set_params({item: params[item]})
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xas_localize_print_program_run_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_localize_print_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_localize_print_program_run_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xas_localize_print_total_dipole_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_localize_print_total_dipole:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_localize_print_total_dipole_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&TOTAL_DIPOLE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END TOTAL_DIPOLE\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xas_localize_print_wannier_centers_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_localize_print_wannier_centers:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_localize_print_wannier_centers_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&WANNIER_CENTERS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END WANNIER_CENTERS\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: parmas[item]})
            else:
                pass


class cp2k_dft_xas_localize_print_wannier_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_localize_print_wannier_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_localize_print_wannier_cubes_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&WANNIER_CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END WANNIER_CUBES\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xas_localize_print_wannier_spreads_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_localize_print_wannier_spreads:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_localize_print_wannier_spreads_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&WANNIER_SPREADS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END WANNIER_SPREADS\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xas_localize_print_wannier_states_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_localize_print_wannier_states_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_localize_print_wannier_states_cubes_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t&END CUBES\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            elif item.splti("-")[6] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xas_localize_print_wannier_states_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_localize_print_wannier_states:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cubes = cp2k_dft_xas_localize_print_wannier_states_cubes()
        self.each = cp2k_dft_xas_localize_print_wannier_states_each()
        
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&WANNIER_STATES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.cubes.status == True:
            self.cubes.to_input(fout)
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END WANIER_STATES\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "CUBES":
                self.cubes.set_params({item: params[item]})
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xas_localize_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.loc_restart = cp2k_dft_xas_localize_print_loc_restart()
        self.molecular_dipoles = cp2k_dft_xas_localize_print_molecular_dipoles()
        self.molecular_states = cp2k_dft_xas_localize_print_molecular_states()
        self.program_run_info = cp2k_dft_xas_localize_print_program_run_info()
        self.total_dipole = cp2k_dft_xas_localize_print_total_dipole()
        self.wannier_centers = cp2k_dft_xas_localize_print_wannier_centers()
        self.wannier_cubes = cp2k_dft_xas_localize_print_wannier_cubes()
        self.wannier_spreads = cp2k_dft_xas_localize_print_wannier_spreads()
        self.wannier_states = cp2k_dft_xas_localize_print_wannier_states()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.loc_restart.status == True:
            self.loc_restart.to_input(fout)
        if self.molecular_dipoles.status == True:
            self.molecular_dipoles.to_input(fout)
        if self.molecular_states.status == True:
            self.molecular_states.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.total_dipole.status == True:
            self.total_dipole.to_input(fout)
        if self.wannier_centers.status == True:
            self.wannier_centers.to_input(fout)
        if self.wannier_cubes.status == True:
            self.wannier_cubes.to_input(fout)
        if self.wannier_spreads.status == True:
            self.wannier_spreads.to_input(fout)
        if self.wannier_states.status == True:
            self.wannier_states.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "LOC_RESTART":
                self.loc_restart.set_params({item: params[item]})
            elif item.split("-")[4] == "MOLECULAR_DIPOLES":
                self.molecular_dipoles.set_params({item: params[item]})
            elif item.split("-")[4] == "MOLECULAR_STATES":
                self.molecular_states.set_params({item: params[item]})
            elif item.split("-")[4] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[4] == "TOTAL_DIPOLE":
                self.total_dipole.set_params({item: params[item]})
            elif item.split("-")[4] == "WANNIER_CENTERS":
                self.wannier_centers.set_params({item: params[item]})
            elif item.split("-")[4] == "WANNIER_CUBES":
                self.wannier_cubes.set_params({item: params[item]})
            elif item.split("-")[4] == "WANNIER_SPREADS":
                self.wannier_spreads.set_params({item: params[item]})
            elif item.split("-")[4] == "WANNIER_STATES":
                self.wannier_states.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xas_localize:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.printout = cp2k_dft_xas_localize_print()
        
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&LOCALIZE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t&END LOCALIZE\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_xas_print_cls_function_cubes_each:
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
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_print_cls_function_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_print_cls_function_cubes_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&CLS_FUNCTION_CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END CLS_FUNCTION_CUBES\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_xas_print_iteration_info_each:
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
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_print_iteration_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_print_iteration_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&ITERATION_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END ITERATION_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_xas_print_loc_restart_each:
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
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_print_loc_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_print_loc_restart_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&LOC_RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
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

class cp2k_dft_xas_print_pdos_each:
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
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_print_pdos_ldos:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&LDOS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END LDOS\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_print_pdos_r_ldos:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&R_LDOS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END R_LDOS\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_print_pdos:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.each = cp2k_dft_xas_print_pdos_each()
        self.ldos = cp2k_dft_xas_print_pdos_ldos()
        self.r_ldos = cp2k_dft_xas_print_pdos_r_ldos()

        # basic seting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PDOS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        if self.ldos.status == True:
            self.ldos.to_input(fout)
        if self.r_ldos.status == True:
            sel.fr_ldos.to_input(fout)  
        fout.write("\t\t\t\t&END PDOS\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            elif item.split("-")[4] == "LDOS":
                self.ldos.set_params({item: params[item]})
            elif item.split("-")[4] == "R_LDOS":
                self.r_ldos.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_xas_print_program_run_info_each:
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
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_print_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.eacj = cp2k_dft_xas_print_program_run_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == True:
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_xas_print_restart_each:
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
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_print_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_print_restart_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END RESTART\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xas_print_wannier_centers_each:
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
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_print_wannier_centers:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_print_wannier_centers_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&WANNIER_CENTERS\n")
        for item in self.params:
            if self.params[item] is not None:
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

class cp2k_dft_xas_print_wannier_cubes_each:
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
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_print_wannier_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_print_wannier_cubes_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&WANNIER_CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
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

class cp2k_dft_xas_print_wannier_spreads_each:
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
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_print_wannier_spreads:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_print_wannier_spreads_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&WANNIER_SPREADS\n")
        for item in self.params:
            if self.params[item] is not None:
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


class cp2k_dft_xas_print_xas_spectrum_each:
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
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_print_xas_spectrum:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_print_xas_spectrum_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&XAS_SPECTRUM\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END XAS_SPECTRUM\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_xas_print_xes_spectrum_each:
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
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_print_xes_spectrum:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_print_xes_spectrum_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&XES_SPECTRUM\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END XES_SPECTRUM\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xas_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cls_function_cubes = cp2k_dft_xas_print_cls_function_cubes()
        self.iteration_info = cp2k_dft_xas_print_iteration_info()
        self.loc_restart = cp2k_dft_xas_print_loc_restart()
        self.pdos = cp2k_dft_xas_print_pdos()
        self.program_run_info = cp2k_dft_xas_print_program_run_info()
        self.restart = cp2k_dft_xas_print_restart()
        self.wanier_centers = cp2k_dft_xas_print_wannier_centers()
        self.wannier_cubes = cp2k_dft_xas_print_wannier_cubes()
        self.wannier_spreads = cp2k_dft_xas_print_wannier_spreads()
        self.xas_spectrum = cp2k_dft_xas_print_xas_spectrum()
        self.xes_spectrum = cp2k_dft_xas_print_xes_spectrum()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.cls_function_cubes.status == True:
            self.cls_function_cubes.to_input(fout)
        if self.iteration_info.status == True:
            self.iteration_info.to_input(fout)
        if self.loc_restart.status == True:
            self.loc_restart.to_input(fout)
        if self.pdos.status == True:
            self.pdos.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.restart.status == True:
            self.restart.to_input(fout)
        if self.wannier_centers.status == True:
            self.wannier_centers.to_input(fout)
        if self.wannier_cubes.status == True:
            self.waniner_cubes.to_input(fout)
        if self.wannier_spreads.status == True:
            self.wannier_spreads.to_input(fout)
        if self.xas_spectrum.status == True:
            self.xas_spectrum.to_input(fout)
        if self.xes_spectrum.status == True:
            self.xes_spectrum.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "CLS_FUNCTION_CUBES":
                self.cls_function_cubes.set_params({item: params[item]})
            elif item.split("-")[3] == "ITERATION_INFO":
                self.iteration_info.set_params({item: params[item]})
            elif item.split("-")[3] == "LOC_RESTART":
                self.loc_restart.set_params({item: params[item]})
            elif item.split("-")[3] == "PDOS":
                self.pdos.set_params({item: params[item]})
            elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[3] == "RESTART":
                self.restart.set_params({item: params[item]})
            elif item.split("-")[3] == "WANNIER_CENTERS":
                self.wannier_centers.set_params({item: params[item]})
            elif item.split("-")[3] == "WANNIER_CUBES":
                self.wannier_cubes.set_params({item: params[item]})
            elif item.split("-")[3] == "WANNIER_SPREADS":
                self.wannier_spreads.set_params({item: params[item]})
            elif item.split("-")[3] == "XAS_SPECTRUM":
                self.xas_spectrum.set_params({item: params[item]})
            elif item.split("-")[3] == "XES_SPECTRUM":
                self.xes_spectrum.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_xas_scf_diagonalization_davidson:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&DAVIDSON\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END DAVIDSON\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_scf_diagonalization_diag_sub_scf_mixing:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&MIXING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END MIXING\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_scf_diagonalization_diag_sub_scf:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.mixing = cp2k_dft_xas_scf_diagonalization_diag_sub_scf_mixing()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&DIAG_SUB_SCF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.mixing.status == True:
            self.mixing.to_input(fout)
        fout.write("\t\t\t\t\t&END DIAG_SUB_SCF\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "MIXING":
                self.mixing.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xas_scf_diagonalization_filter_matrix:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&FILTER_MATRIX\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END FILTER_MATRIX\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_scf_diagonalization_krylov:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&KRYLOV\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END KRYLOV\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_scf_diagonalization_ot:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&OT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END OT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_scf_diagonalization:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.davidson = cp2k_dft_xas_scf_diagonalization_davidson()
        self.diag_sub_scf = cp2k_dft_xas_scf_diagonalization_diag_sub_scf()
        self.filter_matrix = cp2k_dft_xas_scf_diagonalization_filter_matrix()
        self.krylov = cp2k_dft_xas_scf_diagonalization_krylov()
        self.ot = cp2k_dft_xas_scf_diagonalization_ot()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&DIAGONALIZATION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.davidson.status == True:
            self.davidson.to_input(fout)
        if self.diag_sub_scf.status == True:
            self.diag_sub_scf.to_input(fout)
        if self.filter_matrix.status == True:
            self.filter_matrix.to_input(fout)
        if self.krylov.status == True:
            self.krylov.to_input(fout)
        if self.ot.status == True:
            self.ot.to_input(fout)
        fout.write("\t\t\t\t&END DIAGONALIZATION\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "DAVIDSON":
                self.davidson.set_params({item: params[item]})
            elif item.split("-")[4] == "DAIG_SUB_SCF":
                self.diag_sub_scf.set_params({item: params[item]})
            elif item.split("-")[4] == "FILTER_MATRIX":
                self.filter_matrix.set_params({item: params[item]})
            elif item.split("-")[4] == "KRYLOV":
                self.krylov.set_params({item: params[item]})
            elif item.split("-")[4] == "OT":
                self.ot.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xas_scf_mixing:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&MIXING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END MIXING\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_scf_mom:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&MOM\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END MOM\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_scf_ot:
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
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END OT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_scf_outer_scf_cdft_opt:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&CDFT_OPT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END CDFT_OPT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_scf_outer_scf:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cdft_opt = cp2k_dft_xas_scf_outer_scf_cdft_opt()
        
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&OUTER_SCF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.cdft_opt.status == True:
            self.cdft_opt.to_input(fout)
        fout.write("\t\t\t\t&END OUTER_SCF\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CDFT_OPT":
                self.cdft_opt.set_params({item: prams[item]})
            else:
                pass


class cp2k_dft_xas_scf_print_davidson_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_scf_print_davidson:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.each = cp2k_dft_xas_scf_print_davidson_each()
        
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&DAVIDSON\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END DAVIDSON\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({ietm: params[item]})
            else:
                pass

class cp2k_dft_xas_scf_print_detailed_energy_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_scf_print_detailed_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_scf_print_detailed_energy_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&DETAILED_ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END DETAILED_ENERGY\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({ietm: params[item]})
            else:
                pass

class cp2k_dft_xas_scf_print_diag_sub_scf_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_scf_print_diag_sub_scf:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_scf_print_diag_sub_scf_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&DIAG_SUB_SCF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END DIAG_SUB_SCF\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({ietm: params[item]})
            else:
                pass


class cp2k_dft_xas_scf_print_diis_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_scf_print_diis_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each =cp2k_dft_xas_scf_print_diis_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&DIIS_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END DIIS_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({ietm: params[item]})
            else:
                pass

class cp2k_dft_xas_scf_print_filter_matrix_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_scf_print_filter_matrix:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_scf_print_filter_matrix_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&FILTER_MATRIX\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END FILTER_MATRIX\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({ietm: params[item]})
            else:
                pass

class cp2k_dft_xas_scf_print_iteration_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_scf_print_iteration_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_scf_print_iteration_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&ITERATION_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END ITERATION_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({ietm: params[item]})
            else:
                pass


class cp2k_dft_xas_scf_print_lanczos_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_scf_print_lanczos:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_scf_print_lanczos_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&LANCZOS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END LANCZOS\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({ietm: params[item]})
            else:
                pass

class cp2k_dft_xas_scf_print_mos_molden_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_scf_print_mos_molden:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_scf_print_mos_molden_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&MOS_MOLDEN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END MOS_MOLDEN\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({ietm: params[item]})
            else:
                pass

class cp2k_dft_xas_scf_print_mo_magnitude_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_scf_print_mo_magnitude:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_scf_print_mo_magnitude_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&MO_MAGNITUDE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END MO_MAGNITUDE\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({ietm: params[item]})
            else:
                pass


class cp2k_dft_xas_scf_print_mo_orthonormality_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_scf_print_mo_orthonormality:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_scf_print_mo_orthonormality_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&MO_ORTHONORMALITY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END ORTHONORMALITY\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({ietm: params[item]})
            else:
                pass

class cp2k_dft_xas_scf_print_program_run_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_scf_print_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False
        

        self.each = cp2k_dft_xas_scf_print_program_run_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({ietm: params[item]})
            else:
                pass

class cp2k_dft_xas_scf_print_restart_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_scf_print_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_scf_print_restart_each()

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END RESTART\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({ietm: params[item]})
            else:
                pass

class cp2k_dft_xas_scf_print_restart_history_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_scf_print_restart_history:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_scf_print_restart_history_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&RESTART_HISTORY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END RESTART_HISTORY\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({ietm: params[item]})
            else:
                pass

class cp2k_dft_xas_scf_print_total_densities_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xas_scf_print_total_densities:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xas_scf_print_total_densities_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&TOTAL_DENSITIES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END TOTAL_DENSITIES\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({ietm: params[item]})
            else:
                pass

class cp2k_dft_xas_scf_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.davidson = cp2k_dft_xas_scf_print_davidson()
        self.detailed_energy = cp2k_dft_xas_scf_print_detailed_energy()
        self.diag_sub_scf = cp2k_dft_xas_scf_print_diag_sub_scf()
        self.diis_info = cp2k_dft_xas_scf_print_diis_info()
        self.filter_matrix = cp2k_dft_xas_scf_print_filter_matrix()
        self.iteration_info = cp2k_dft_xas_scf_print_iteration_info()
        self.lanczos = cp2k_dft_xas_scf_print_lanczos()
        self.mos_molden = cp2k_dft_xas_scf_print_mos_molden()
        self.mo_magnitude = cp2k_dft_xas_scf_print_mo_magnitude()
        self.mo_orthonormality = cp2k_dft_xas_scf_print_mo_orthonormality()
        self.program_run_info = cp2k_dft_xas_scf_print_program_run_info()
        self.restart = cp2k_dft_xas_scf_print_restart()
        self.restart_history = cp2k_dft_xas_scf_print_restart_history()
        self.total_densities = cp2k_dft_xas_scf_print_total_densities()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
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
        if self.iteration_info.status == True:
            self.iteration_info.to_input(fout)
        if self.lanczos.status == True:
            self.lanczos.to_input(fout)
        if self.mos_molden.status == True:
            self.mos_molden.to_input(fout)
        if self.mo_magnitude.status == True:
            self.mo_magnitude.to_input(fout)
        if self.mo_orthonormality.status == True:
            self.mo_orthonormality.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.restart.status == True:
            self.restart.to_input(fout)
        if self.restart_history.status == True:
            self.restart_history.to_input(fout)
        if self.total_densities.status == True:
            self.total_densities.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "DAVIDSON":
                self.davidson.set_params({item: params[item]})
            elif item.split("-")[4] == "DETAILED_ENERGY":
                self.detailed_energy.set_params({item: params[item]})
            elif item.split("-")[4] == "DIAG_SUB_SCF":
                self.diag_sub_scf.set_params({item: params[item]})
            elif item.split("-")[4] == "DIIS_INFO":
                self.diis_info.set_params({item: params[item]})
            elif item.split("-")[4] == "FILTER_MATRIX":
                self.filter_matrix.set_params({item: params[item]})
            elif item.split("-")[4] == "ITERATION_INFO":
                self.iteration_info.set_params({item: params[item]})
            elif item.split("-")[4] == "LANCZOS":
                self.lanczos.set_params({item: params[item]})
            elif item.split("-")[4] == "MOS_MOLDEN":
                self.mos_molden.set_params({item: params[item]})
            elif item.split("-")[4] == "MO_MAGNITUDE":
                self.mo_magnitude.set_params({item: params[item]})
            elif item.split("-")[4] == "MO_ORTHONORMALITY":
                self.mo_orthonormality.set_params({item: params[item]})
            elif item.split("-")[4] == "PROGRAM_RUN_INFO":
                self.proram_run_info.set_params({item: params[item]})
            elif item.split("-")[4] == "RESTART":
                self.restart.set_params({item: params[item]})
            elif item.split("-")[4] == "RESTART_HISTORY":
                self.restart_history.set_params({item: params[item]})
            elif item.split("-")[4] == "TOTAL_DENSITIES":
                self.total_densities.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xas_scf_smear:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&SMEAR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END SMEAR\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xas_scf:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.diagonalization = cp2k_dft_xas_scf_diagonalization()
        self.mixing = cp2k_dft_xas_scf_mixing()
        self.mom = cp2k_dft_xas_scf_mom()
        self.ot = cp2k_dft_xas_scf_ot()
        self.outer_scf = cp2k_dft_xas_scf_outer_scf()
        self.printout = cp2k_dft_xas_scf_print()
        self.smear = cp2k_dft_xas_scf_smear()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&SCF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.diagonalization.status == True:
            self.diagonalization.to_input(fout)
        if self.mixing.status == True:
            self.mixing.to_input(fout)
        if self.mom.status == True:
            self.mom.to_input(fout)
        if self.ot.status == True:
            self.ot.to_input(fout)
        if self.outer_scf.status == True:
            self.outer_scf.to_input(fout)
        if self.printout.statu == True:
            self.printout.to_input(fout)
        if self.smear.status == True:
            self.smear.to_input(fout)
        fout.write("\t\t\t&END SCF\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "DIAGONALIZATION":
                self.diagonalization.set_params({item: params[item]})
            elif item.split("-")[3] == "MIXING":
                self.mixing.set_params({item: params[item]})
            elif item.split("-")[3] == "MOM":
                self.mom.set_params({item: params[item]})
            elif item.split("-")[3] == "OT":
                self.ot.set_params({item: params[item]})
            elif item.split("-")[3] == "OUTER_SCF":
                self.outer_scf.set_params({item: params[item]})
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[3] == "SMEAR":
                self.smear.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_xas:
    def __init__(self):
        self.params = {
                "ADDED_MOS": None,
                "METHOD": None,
                }
        self.status = False

        self.localize = cp2k_dft_xas_localize()
        self.printout = cp2k_dft_xas_print()
        self.scf = cp2k_dft_xas_scf()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&XAS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.localize.status == True:
            self.localize.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.scf.status == True:
            self.scf.to_input(fout)
        fout.write("\t\t&END XAS\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "LOCALIZE":
                self.localize.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[2] == "SCF":
                self.scf.set_params({item: params[item]})
            else:
                pass
