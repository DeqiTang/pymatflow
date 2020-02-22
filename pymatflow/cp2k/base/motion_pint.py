#!/usr/bin/env python
# _*_ coding: utf-8 _*_


class cp2k_motion_pint_beads_coord:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&COORD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END COORD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_pint_beads_velocity:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&VELOCITY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END VELOCITY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_beads:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.coord = cp2k_motion_pint_beads_coord()
        self.velocity = cp2k_motion_pint_beads_velocity()


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&BEADS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.coord.status == True:
            self.coord.to_input(fout)
        if self.velocity.status == True:
            self.velocity.to_input(fout)
        fout.write("\t\t&END BEADS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "COORD":
                self.coord.set_params({item: params[item]})
            elif item.split("-")[2] == "VELOCITY":
                self.velocity.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_gle_rng_init:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RNG_INIT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RNG_INIT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_gle_s:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&S\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END S\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_pint_gle_thermostat_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&THERMOSTAT_ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END THERMOSTAT_ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_gle:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.rng_init = cp2k_motion_pint_gle_rng_init()
        self.s = cp2k_motion_pint_gle_s()
        self.thermostat_energy = cp2k_motion_pint_gle_thermostat_energy()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&GLE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.rng_init.status == True:
            self.rng_init.to_input(fout)
        if self.s.status == True:
            self.s.to_input(fout)
        if self.thermostat_energy.status == True:
            self.thermostat_energy.to_input(fout)
        fout.write("\t\t&END GLE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "RNG_INIT":
                self.rng_init.set_params({item: params[item]})
            elif item.split("-")[2] == "S":
                self.s.set_params({item: params[item]})
            elif item.split("-")[2] == "THERMOSTAT_ENERGY":
                self.thermostat_energy.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_pint_helium_averages:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&AVERAGES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END AVERAGES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_ceperley_m_sampling:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&M-SAMPLING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END M-SAMPLING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_ceperley:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.m_sampling = cp2k_motion_pint_helium_ceperley_m_sampling()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CEPERLEY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.m_sampling.status == True:
            self.m_sampling.to_input(fout)
        fout.write("\t\t\t&END CEPERLEY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "M_SAMPLING":
                self.m_sampling.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_pint_helium_coord:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&COORD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END COORD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_force:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&FORCE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END FORCE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_perm:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PERM\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END PERM\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_motion_pint_helium_print_accepts_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_accepts:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_accepts_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&ACCEPTS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END ACCEPTS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_helium_print_action_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_action:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_action_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&ACTION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END ACTION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_pint_helium_print_coordinates_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_coordinates:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_coordinates_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&COORDINATES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END COORDINATES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_pint_helium_print_energy_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_energy_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_helium_print_forces_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_forces:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_forces_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&FORCES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END FORCES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_helium_print_forces_inst_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_forces_inst:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_forces_inst_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&FORCES_INST\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END FORCES_INST\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_helium_print_moment_of_inertia_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_moment_of_inertia:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_moment_of_inertia_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&MOMENT_OF_INERTIA\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END MOMENT_OF_INERTIA\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_helium_print_moment_of_inertia_avg_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_moment_of_inertia_avg:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_moment_of_inertia_avg_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&MOMENT_OF_INERTIA_AVG\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END MOMENT_OF_INERTIA_AVG\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_pint_helium_print_perm_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_perm:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_perm_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PERM\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END PERM\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_pint_helium_print_plength_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_plength:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_plength_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PLENGTH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END PLENGTH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_pint_helium_print_projected_area_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_projected_area:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_projected_area_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PROJECTED_AREA\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END PROJECTED_AREA\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_pint_helium_print_projected_area_2_avg_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_projected_area_2_avg:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_projected_area_2_avg_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PROJECTED_AREA_2_AVG\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END PROJECTED_AREA_2_AVG\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_pint_helium_print_rdf_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_rdf:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_rdf_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&RDF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END RDF\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_pint_helium_print_rho_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_rho:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_rho_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&RHO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END RHO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_helium_print_winding_number_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_winding_number:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_winding_number_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&WINDING_NUMBER\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END WINDING_NUMBER\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_pint_helium_print_winding_number_2_avg_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_print_winding_number_2_avg:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_helium_print_winding_number_2_avg_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&WINDING_NUMBER_2_AVG\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END WINDING_NUMBER_2_AVG\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_pint_helium_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.accepts = cp2k_motion_pint_helium_print_accepts()
        self.action = cp2k_motion_pint_helium_print_action()
        self.coordinates = cp2k_motion_pint_helium_print_coordinates()
        self.energy = cp2k_motion_pint_helium_print_energy()
        self.forces = cp2k_motion_pint_helium_print_forces()
        self.forces_inst = cp2k_motion_pint_helium_print_forces_inst()
        self.moment_of_inertia = cp2k_motion_pint_helium_print_moment_of_inertia()
        self.moment_of_inertia_avg = cp2k_motion_pint_helium_print_moment_of_inertia_avg()
        self.perm = cp2k_motion_pint_helium_print_perm()
        self.plength = cp2k_motion_pint_helium_print_plength()
        self.projected_area = cp2k_motion_pint_helium_print_projected_area()
        self.projected_area_2_avg = cp2k_motion_pint_helium_print_projected_area_2_avg()
        self.rdf = cp2k_motion_pint_helium_print_rdf()
        self.rho = cp2k_motion_pint_helium_print_rho()
        self.winding_number = cp2k_motion_pint_helium_print_winding_number()
        self.winding_number_2_avg = cp2k_motion_pint_helium_print_winding_number_2_avg()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.accepts.status == True:
            self.accepts.to_input(fout)
        if self.action.status == True:
            self.action.to_input(fout)
        if self.coordinates.status == True: 
            self.coordinates.to_input(fout)
        if self.energy.status == True: 
            self.energy.to_input(fout)
        if self.forces.status == True:
            self.forces.to_input(fout)
        if self.forces_inst.status == True: 
            self.forces_inst.to_input(fout)
        if self.moment_of_inertia.status == True: 
            self.moment_of_inertia.to_input(fout)
        if self.moment_of_inertia_avg.status == True: 
            self.moment_of_inertia_avg.to_input(fout)
        if self.perm.status == True: 
            self.perm.to_input(fout)
        if self.plength.status == True: 
            self.plength.to_input(fout)
        if self.projected_area.status == True: 
            self.projected_area.to_input(fout)
        if self.projected_area_2_avg.status == True: 
            self.projected_area_2_avg.to_input(fout)
        if self.rdf.status == True:
            self.rdf.to_input(fout)
        if self.rho.status == True:
            self.rho.to_input(fout)
        if self.winding_number.status == True: 
            self.winding_number.to_input(fout)
        if self.winding_number_2_avg.status == True: 
            self.winding_number_2_avg.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "ACCEPTS":
                self.accepts.set_params({item: params[item]})
            elif item.split("-")[3] == "ACTION":
                self.action.to_input(fout)
            elif item.split("-")[3] == "COORDINATES":
                self.coordinates.to_input(fout)
            elif item.split("-")[3] == "ENERGY":
                self.energy.to_input(fout)
            elif item.split("-")[3] == "FORCES":
                self.forces.to_input(fout)
            elif item.split("-")[3] == "FORCES_INST":
                self.forces_inst.to_input(fout)
            elif item.split("-")[3] == "MOMENT_OF_INERTIA":
                self.moment_of_inertia.to_input(fout)
            elif item.split("-")[3] == "MOMENT_OF_INERTIA_AVG":
                self.moment_of_inertia_avg.to_input(fout)
            elif item.split("-")[3] == "PERM":
                self.perm.to_input(fout)
            elif item.split("-")[3] == "PLENGTH":
                self.plength.to_input(fout)
            elif item.split("-")[3] == "PROJECTED_AREA":
                self.projected_area.to_input(fout)
            elif item.split("-")[3] == "PROJECTED_AREA_2_AVG":
                self.projected_area_2_avg.to_input(fout)
            elif item.split("-")[3] == "RDF":
                self.rdf.to_input(fout)
            elif item.split("-")[3] == "RHO":
                self.rho.to_input(fout)
            elif item.split("-")[3] == "WINDING_NUMBER":
                self.winding_number.to_input(fout)
            elif item.split("-")[3] == "WINDING_NUMBER_2_AVG":
                self.winding_number_2_avg.to_input(fout)
            else:
                pass


class cp2k_motion_pint_helium_rdf:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RDF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RDF\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_pint_helium_rho:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RHO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RHO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_rng_state:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RNG_STATE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RNG_STATE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium_worm:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&WORM\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END WORM\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_helium:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.avearges = cp2k_motion_pint_helium_averages()
        self.ceperley = cp2k_motion_pint_helium_ceperley()
        self.coord = cp2k_motion_pint_helium_coord()
        self.force = cp2k_motion_pint_helium_force()
        self.perm = cp2k_motion_pint_helium_perm()
        self.printout = cp2k_motion_pint_helium_print()
        self.rdf = cp2k_motion_pint_helium_rdf()
        self.rho = cp2k_motion_pint_helium_rho()
        self.rng_state = cp2k_motion_pint_helium_rng_state()
        self.worm = cp2k_motion_pint_helium_worm()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&HELIUM\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.averages.status == True:
            self.averages.to_input(fout)
        if self.ceperley.status == True:
            self.ceperley.to_input(fout)
        if self.coord.status == True:
            self.coord.to_input(fout)
        if self.force.status == True:
            self.force.to_input(fout)
        if self.perm.status == True:
            self.perm.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.rdf.status == True:
            self.rdf.to_input(fout)
        if self.rho.status == True:
            self.rho.to_input(fout)
        if self.rng_state.status == True:
            self.rng_state.to_input(fout)
        if self.worm.status == True:
            self.worm.to_input(fout)
        fout.write("\t\t&END HELIUM\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "AVERAGES":
                self.averages.set_params({item: params[item]})
            elif item.split("-")[2] == "CEPERLEY":
                self.ceperley.set_params({item: params[item]})
            elif item.split("-")[2] == "COORD":
                self.coord.set_params({item: params[item]})
            elif item.split("-")[2] == "FORCE":
                self.force.set_params({item: params[item]})
            elif item.split("-")[2] == "PERM":
                self.perm.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[2] == "RDF":
                self.rdf.set_params({item: params[item]})
            elif item.split("-")[2] == "RHO":
                self.rho.set_params({item: params[item]})
            elif item.split("-")[2] == "RNG_STATE":
                self.rng_state.set_params({item: params[item]})
            elif item.split("-")[2] == "WORM":
                self.worm.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_pint_init:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&INIT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END INIT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_pint_normalmode:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&NORMALMODE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END NORMALMODE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_pint_nose_coord:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&COORD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END COORD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_pint_nose_velocity:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&VELOCITY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END VELOCITY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_nose:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.coord = cp2k_motion_pint_nose_coord()
        self.velocity = cp2k_motion_pint_nose_velocity()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&NOSE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.coord.status == True:
            self.coord.to_input(fout)
        if self.velocity.status == True:
            self.velocity.to_input(fout)
        fout.write("\t\t&END NOSE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "COORD":
                self.coord.set_params({item: params[item]})
            elif item.split("-")[2] == "VELOCITY":
                self.velocity.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_piglet_extra_dof:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&EXTRA_DOF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END EXTRA_DOF\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_pint_piglet_rng_init:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&RNG_INIT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END RNG_INIT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_piglet:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.extra_dof = cp2k_motion_pint_piglet_extra_dof()
        self.rng_init = cp2k_motion_pint_piglet_rng_init()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PIGLET\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.extra_dof.status == True:
            self.extra_dfo.to_input(fout)
        if self.rng_init.status == True:
            self.rng_init.to_input(fout)
        fout.write("\t\t&END PIGLET\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EXTRA_DOF":
                self.extra_dof.set_params({item: params[item]})
            elif item.split("-")[2] == "RNG_INIT":
                self.rng_init.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_pile_rng_init:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&RNG_INIT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END RNG_INIT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_pint_pile:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.rng_init = cp2k_motion_pint_pile_rng_init()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PILE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.rng_init.status == True:
            self.rng_init.to_input(fout)
        fout.write("\t\t&END PILE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "RNG_INIT":
                self.rng_init.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_pint_print_action_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
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

class cp2k_motion_pint_print_action:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_print_action_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&ACTION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END ACTION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_print_centroid_gyr_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
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

class cp2k_motion_pint_print_centroid_gyr:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_print_centroid_gyr_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CENTROID_GYR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END CENTROID_GYR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_print_centroid_pos_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
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

class cp2k_motion_pint_print_centroid_pos:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_print_centroid_pos_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CENTROID_POS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END CENTROID_POS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_print_centroid_vel_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
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

class cp2k_motion_pint_print_centroid_vel:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_print_centroid_vel_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CENTROID_VEL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END CENTROID_VEL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_print_com_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
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

class cp2k_motion_pint_print_com:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_print_com_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&COM\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END COM\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_print_energy_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
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

class cp2k_motion_pint_print_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_pint_print_energy_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.action = cp2k_motion_pint_print_action()
        self.centroid_gyr = cp2k_motion_pint_print_centroid_gyr()
        self.centroid_pos = cp2k_motion_pint_print_centroid_pos()
        self.centroid_vel = cp2k_motion_pint_print_centroid_vel()
        self.com = cp2k_motion_pint_print_com()
        self.energy = cp2k_motion_pint_print_energy()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.action.status == True:
            slef.action.to_input(fout)
        if self.centroid_gyr.status == True:
            self.centroid_gyr.to_input(fout)
        if slef.centroid_pos.status == True:
            self.centroid_pos.to_input(fout)
        if self.centroid_vel.status == True:
            self.centroid_vel.to_input(fout)
        if self.com.status == True:
            self.com.to_input(fout)
        if self.energy.status == True:
            self.energy.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "ACTION":
                self.action.set_params({item: params[item]})
            elif item.split("-")[2] == "CENTROID_GYR":
                self.centroid_gyr.set_params({item: params[item]})
            elif item.split("-")[2] == "CENTROID_POS":
                self.centroid_pos.set_params({item: params[item]})
            elif item.split("-")[2] == "CENTROID_VEL":
                self.centroid_vel.set_params({item: params[item]})
            elif item.split("-")[2] == "COM":
                self.com.set_params({item: params[item]})
            elif item.split("-")[2] == "ENERGY":
                self.energy.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_pint_qtb_rng_init:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RNG_INIT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RNG_INIT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_pint_qtb:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.rng_init = cp2k_motion_pint_qtb_rng_init()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&QTB\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.rng_init.status == True:
            self.rng_init.to_input(fout)
        fout.write("\t\t&END QTB\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "RNG_INIT":
                self.rng_init.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_pint_staging:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&STAGING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END STAGING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_motion_pint:
    def __init__(self):
        self.params = {
                "DT": None,
                "FIX_CENTROID_POS": None,
                "HARM_INT": None,
                "ITERATION": None,
                "MAX_STEP": None,
                "NRESPA": None,
                "NUM_STEPS": None,
                "P": None,
                "PROC_PER_REPLICA": None,
                "PROPAGATOR": None,
                "TEMP": None,
                "TRANSFORMATION": None,
                "T_TOL": None,
                }
        self.status = False

        self.beads = cp2k_motion_pint_beads()
        self.gle = cp2k_motion_pint_gle()
        self.helium = cp2k_motion_pint_helium()
        self.init = cp2k_motion_pint_init()
        self.normalmode = cp2k_motion_pint_normalmode()
        self.nose = cp2k_motion_pint_nose()
        self.piglet = cp2k_motion_pint_piglet()
        self.pile = cp2k_motion_pint_pile()
        self.printout = cp2k_motion_pint_print()
        self.qtb = cp2k_motion_pint_qtb()
        self.staging = cp2k_motion_pint_staging()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&PINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.beads.status == True:
            self.beads.to_input(fout)
        if self.gle.status == True:
            self.gle.to_input(fout)
        if self.helium.status == True:
            self.helium.to_input(fout)
        if self.init.status == True:
            self.init.to_input(fout)
        if self.normalmode.status == True:
            self.normalmode.to_input(fout)
        if self.nose.status == True:
            self.nose.to_input(fout)
        if self.piglet.status == True:
            self.piglet.to_input(fout)
        if self.pile.status == True:
            self.pile.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.qtb.status == True:
            self.qtb.to_input(fout)
        if self.staging.status == True:
            self.staging.to_input(fout)
        fout.write("\t&END PINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-") == "BEADS":
                self.beads.set_params({item: params[item]})
            elif item.split("-") == "GLE":
                self.gle.set_params({item: params[item]})
            elif item.split("-") == "HELIUM":
                self.helium.set_params({item: params[item]})
            elif item.split("-") == "INIT":
                self.init.set_params({item: params[item]})
            elif item.split("-") == "NORMALMODE":
                self.normalmode.set_params({item: params[item]})
            elif item.split("-") == "NOSE":
                self.nose.set_params({item: params[item]})
            elif item.split("-") == "PIGLET":
                self.piglet.set_params({item: params[item]})
            elif item.split("-") == "PILE":
                self.pile.set_params({item: params[item]})
            elif item.split("-") == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-") == "QTB":
                self.qtb.set_params({item: params[item]})
            elif item.split("-") == "STAGING":
                self.staging.set_params({item: params[item]})
            else:
                pass

