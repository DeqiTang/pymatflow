#!/usr/bin/env python
# _*_ coding: utf-8 _*_


class cp2k_motion_free_energy_alchemical_change:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&ALCHEMICAL_CHANGE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END ALCHEMICAL_CHANGE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_free_energy_free_energy_info_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_free_energy_free_energy_info:
    def __init__(self):
        self.params = {}
        self.status = False

        self.each = cp2k_motion_free_energy_free_energy_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&FREE_ENERGY_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END FREE_ENERGY_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_free_energy_metadyn_ext_lagrange_fs:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EXT_LAGRANGE_FS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END EXT_LAGRANGE_FS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_free_energy_metadyn_ext_lagrange_ss:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EXT_LAGRANGE_SS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END EXT_LAGRANGE_SS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_motion_free_energy_metadyn_ext_lagrange_ss0:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EXT_LAGRANGE_SS0\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END EXT_LAGRANGE_SS0\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_free_energy_metadyn_ext_lagrange_vvp:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EXT_LAGRANGE_VVP\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END EXT_LAGRANGE_VVP\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_free_energy_metadyn_metavar_wall_gaussian:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&GAUSSIAN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END GAUSSIAN\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_free_energy_metadyn_metavar_wall_quadratic:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&QUADRATIC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END QUADRATIC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_free_energy_metadyn_metavar_wall_quartic:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&QUARTIC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END QUARTIC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_free_energy_metadyn_metavar_wall_reflective:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&REFLECTIVE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END REFLECTIVE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_free_energy_metadyn_metavar_wall:
    def __init__(self):
        self.params = {}
        self.status = False

        self.gaussian = cp2k_motion_free_energy_metadyn_metavar_wall_gaussian()
        self.quadratic = cp2k_motion_free_energy_metadyn_metavar_wall_quadratic()
        self.quartic = cp2k_motion_free_energy_metadyn_metavar_wall_quartic()
        self.reflective = cp2k_motion_free_energy_metadyn_metavar_wall_reflective()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&WALL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.gaussian.status == True:
            self.gaussian.to_input(fout) 
        if self.quadratic.status == True:
            self.quadratic.to_input(fout)
        if self.quartic.status == True:
            self.quartic.to_input(fout)
        if self.reflective.status == True:
            self.reflective.to_input(fout)
        fout.write("\t\t\t\t&END WALL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "GAUSSIAN":
                self.gaussian.set_params({item: params[item]})
            elif item.split("-")[4] == "QUADRATIC":
                self.quadratic.set_params({item: params[item]})
            elif item.split("-")[4] == "QUARTIC":
                self.quartic.set_params({item: params[item]})
            elif item.split("-")[4] == "REFLECTIVE":
                self.reflective.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_free_energy_metadyn_metavar:
    def __init__(self):
        self.params = {}
        self.status = False

        self.wall = cp2k_motion_free_energy_metadyn_metavar_wall()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&METAVAR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.wall.status == True:
            self.wall.to_input(fout)
        fout.write("\t\t\t&END METAVAR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "WALL":
                self.wall.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_free_energy_metadyn_multiple_walkers_walkers_file_name:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&WALKERS_FILE_NAME\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END WALKERS_FILE_NAME\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_free_energy_metadyn_multiple_walkers:
    def __init__(self):
        self.params = {}
        self.status = False

        self.walkers_file_name = cp2k_motion_free_energy_metadyn_multiple_walkers_walkers_file_name()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&MULTIPLE_WALKERS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.walkers_file_name.status == True:
            self.walkers_file_name.to_input(fout)
        fout.write("\t\t\t&END MULTIPLE_WALKERS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "WALKERS_FILE_NAME":
                self.walkers.file_name.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_free_energy_metadyn_print_colvar_each:
    def __init__(self):
        self.params = {}
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
            else:
                pass


class cp2k_motion_free_energy_metadyn_print_colvar:
    def __init__(self):
        self.params = {}
        self.status = False

        self.each = cp2k_motion_free_energy_metadyn_print_colvar_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&COLVAR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END COLVAR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_free_energy_metadyn_print_hills_each:
    def __init__(self):
        self.params = {}
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
            else:
                pass


class cp2k_motion_free_energy_metadyn_print_hills:
    def __init__(self):
        self.params = {}
        self.status = False

        self.each = cp2k_motion_free_energy_metadyn_print_hills_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&HILLS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END HILLS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_free_energy_metadyn_print_program_run_info_each:
    def __init__(self):
        self.params = {}
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
            else:
                pass


class cp2k_motion_free_energy_metadyn_print_program_run_info:
    def __init__(self):
        self.params = {}
        self.status = False

        self.each = cp2k_motion_free_energy_metadyn_print_program_run_info_each()
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
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_free_energy_metadyn_print_temperature_colvar_each:
    def __init__(self):
        self.params = {}
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
            else:
                pass


class cp2k_motion_free_energy_metadyn_print_temperature_colvar:
    def __init__(self):
        self.params = {}
        self.status = False

        self.each = cp2k_motion_free_energy_metadyn_print_temperature_colvar_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&TEMPERATURE_COLVAR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END TEMPERATURE_COLVAR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_free_energy_metadyn_print:
    def __init__(self):
        self.params = {}
        self.status = False

        self.colvar = cp2k_motion_free_energy_metadyn_print_colvar()
        self.hills = cp2k_motion_free_energy_metadyn_print_hills()
        self.program_run_info = cp2k_motion_free_energy_metadyn_print_program_run_info()
        self.temperature_colvar = cp2k_motion_free_energy_metadyn_print_temperature_colvar()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.colvar.status == True:
            self.colvar.to_input(fout)
        if self.hills.status == True:
            self.hills.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.temperature_colvar.status == True:
            self.temperature_colvar.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "COLVAR":
                self.colvar.set_params({item: params[item]})
            elif item.split("-")[3] == "HILLS":
                self.hills.set_params({item: params[item]})
            elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[3] == "TEMPERATURE_COLVAR":
                self.temperature_colvar.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_free_energy_metadyn_spawned_hills_height:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&SPAWNED_HILLS_HEIGHT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END SPAWNED_HILLS_HEIGHT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_free_energy_metadyn_spawned_hills_invdt:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&SPAWNED_HILLS_INVDT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END SPAWNED_HILLS_INVDT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_free_energy_metadyn_spawned_hills_pos:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&SPAWNED_HILLS_POS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END SPAWNED_HILLS_POS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_free_energy_metadyn_spawned_hills_scale:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&SPAWNED_HILLS_SCALE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END SPAWNED_HILLS_SCALE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_free_energy_metadyn:
    def __init__(self):
        self.params = {}
        self.status = False


        self.ext_lagrange_fs = cp2k_motion_free_energy_metadyn_ext_lagrange_fs()
        self.ext_lagrange_ss = cp2k_motion_free_energy_metadyn_ext_lagrange_ss()
        self.ext_lagrange_ss0 = cp2k_motion_free_energy_metadyn_ext_lagrange_ss0()
        self.ext_lagrange_vvp = cp2k_motion_free_energy_metadyn_ext_lagrange_vvp()
        self.metavar = cp2k_motion_free_energy_metadyn_metavar()
        self.multiple_walkers = cp2k_motion_free_energy_metadyn_multiple_walkers()
        self.printout = cp2k_motion_free_energy_metadyn_print()
        self.spawned_hills_height = cp2k_motion_free_energy_metadyn_spawned_hills_height()
        self.spawned_hills_invdt = cp2k_motion_free_energy_metadyn_spawned_hills_invdt()
        self.spawned_hills_pos = cp2k_motion_free_energy_metadyn_spawned_hills_pos()
        self.spawned_hills_scale = cp2k_motion_free_energy_metadyn_spawned_hills_scale()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&METADYN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.ext_lagrange_fs.status == True:
            self.ext_lagrange_fs.to_input(fout)
        if self.ext_lagrange_ss.status == True:
            self.ext_lagrange_ss.to_input(fout)
        if self.ext_lagrange_ss0.status == True:
            self.ext_lagrange_ss0.to_input(fout)
        if self.ext_lagrange_vvp.status == True:
            self.ext_lagrange_vvp.to_input(fout)
        if self.metavar.status == True:
            self.metavar.to_input(fout)
        if self.multiple_walkers.status == True:
            self.multiple_walkers.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.spawned_hills_height.status == True:
            self.spawned_hills_height.to_input(fout)
        if self.spawned_hills_invdt.status == True:
            self.spawned_hills_invdt.to_input(fout)
        if self.spawned_hills_pos.status == True:
            self.spanwed_hills_pos.to_input(fout)
        if self.spawned_hills_scale.status == True:
            self.spawned_hills_scale.to_input(fout)
        fout.write("\t\t&END METADYN\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EXT_LAGRANGE_FS":
                self.ext_lagrange_fs.set_params({item: params[item]})
            elif item.split("-")[2] == "EXT_LAGRANGE_SS":
                self.ext_lagrange_ss.set_params({item: params[item]})
            elif item.split("-")[2] == "EXT_LAGRANGE_SS0":
                self.ext_lagrange_ss0.set_params({item: params[item]})
            elif item.split("-")[2] == "EXT_LAGRANGE_VVP":
                self.ext_lagrange_vvp.set_params({item: params[item]})
            elif item.split("-")[2] == "METAVAR":
                self.metavar.set_params({item: params[item]})
            elif item.split("-")[2] == "MULTIPLE_WALKERS":
                self.multiple_walkers.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[2] == "SPAWNED_HILLS_HEIGHT":
                self.spawned_hills_height.set_params({item: params[item]})
            elif item.split("-")[2] == "SPAWNED_HILLS_INVDT":
                self.spawned_hills_invdt.set_params({item: params[item]})
            elif item.split("-")[2] == "SPAWNED_HILLS_POS":
                self.spawned_hills_pos.set_params({item: params[item]})
            elif item.split("-")[2] == "SPANWED_HILLS_SCALE":
                self.spawned_hills_scale.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_free_energy_umbrella_integration_convergence_control:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CONVERGENCE_CONTROL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END CONVERGENCE_CONTROL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_free_energy_umbrella_integration_uvar:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&UVAR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END UVAR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_free_energy_umbrella_integration:
    def __init__(self):
        self.params = {}
        self.status = False

        self.convergence_control = cp2k_motion_free_energy_umbrella_integration_convergence_control()
        self.uvar = cp2k_motion_free_energy_umbrella_integration_uvar()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&UMBRELLA_INTEGRATION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.convergence_control.status == True:
            self.convergence_control.to_input(fout)
        if self.uvar.status == True:
            self.uvar.to_input(fout)
        fout.write("\t\t&END UMBRELLA_INTERGRATION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "CONVERGENCE_CONTROL":
                self.convergence_control.set_params({item: params[item]})
            elif item.split("-")[2] == "UVAR":
                self.uvar.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_free_energy:
    def __init__(self):
        self.params = {}
        self.status = False

        self.alchemical_change = cp2k_motion_free_energy_alchemical_change()
        self.free_energy_info = cp2k_motion_free_energy_free_energy_info()
        self.metadyn = cp2k_motion_free_energy_metadyn()
        self.umbrella_integration = cp2k_motion_free_energy_umbrella_integration()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&FREE_ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, str(self.params[item])))
        if self.alchemical_change.status == True:
            self.alchmical_change.to_input(fout)
        if self.free_energy_info.status == True:
            self.free_energy_Info.to_input(fout)
        if self.metadyn.status == True:
            self.metadyn.to_input(fout)
        if self.umbrella_integration.status == True:
            self.umbrella_integration.to_input(fout)
        fout.write("\t&END FREE_ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "ALCHEMICAL_CHANGE":
                self.alchemical_change.set_params({item: params[item]})
            elif item.split("-")[1] == "FREE_ENERGY_INFO":
                slef.free_energy_info.set_params({item: params[item]})
            elif item.split("-")[1] == "METADYN":
                self.metadyn.set_params({item: params[item]})
            elif item.split("-")[1] == "UMBRELLA_INTERGRATION":
                self.umbrella_integration.set_params({item: params[item]})
            else:
                pass

