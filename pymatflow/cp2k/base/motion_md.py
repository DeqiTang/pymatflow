#!/usr/bin/env python
# _*_ coding: utf-8 _*_


class cp2k_motion_md_adiabatic_dynamics_thermostat_fast_define_region:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t&DEFINE_REGION\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END_DEFINE_REGION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_adiabatic_dynamics_thermostat_fast_nose_coord:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&COORD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END COORD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_adiabatic_dynamics_thermostat_fast_nose_force:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&FORCE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END FORCE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_adiabatic_dynamics_thermostat_fast_nose_mass:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&MASS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END MASS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_adiabatic_dynamics_thermostat_fast_nose_velocity:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&VELOCITY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END VELOCITY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_adiabatic_dynamics_thermostat_fast_nose:
    def __init__(self):
        self.params = {
                "LENGTH": None,
                "MTS": None,
                "TIMECON": None,
                "YOSHIDA": None,
                }
        self.status = False

        self.coord = cp2k_motion_md_adiabatic_dynamics_thermostat_fast_nose_coord()
        self.force = cp2k_motion_md_adiabatic_dynamics_thermostat_fast_nose_force()
        self.mass = cp2k_motion_md_adiabatic_dynamics_thermostat_fast_nose_mass()
        self.velocity = cp2k_motion_md_adiabatic_dynamics_thermostat_fast_nose_velocity()
        # basic setting
        self.params["LENGTH"] = 3
        self.params["MTS"] = 2
        self.params["TIMECON"] = 1.0e3
        self.params["YOSHIDA"] = 3

    def to_input(self, fout):
        fout.write("\t\t\t\t&NOSE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.coord.status == True:
            self.coord.to_input(fout)
        if self.force.status == True:
            self.force.to_input(fout)
        if self.mass.status == True:
            self.mass.to_input(fout)
        if self.velocity.status == True:
            self.velocity.to_input(fout)
        fout.write("\t\t\t\t&END NOSE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "VELOCITY":
                self.velocity.set_params({item: params[item]})
            elif item.split("-")[4] == "MASS":
                self.mass.set_params({item: params[item]})
            elif item.split("-")[4] == "FORCE":
                self.force.set_params({item: params[item]})
            elif item.split("-")[4] == "COORD":
                self.coord.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_adiabatic_dynamics_thermostat_fast:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.define_region = cp2k_motion_md_adiabatic_dynamics_thermostat_fast_define_region()
        self.nose = cp2k_motion_md_adiabatic_dynamics_thermostat_fast_nose()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&THERMOSTAT_FAST\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.define_region.status == True:
            self.defin_region.to_input(fout)
        if self.nose.status == True:
            self.nose.to_input(fout)
        fout.write("\t\t\t&END THERMOSTAT_FAST\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "DEFINE_REGION":
                self.define_region.set_params({item: params[item]})
            elif item.split("-")[3] == "NOSE":
                self.nose.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_md_adiabatic_dynamics_thermostat_slow_define_region:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t&define_region\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&end define_region\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_adiabatic_dynamics_thermostat_slow_nose_coord:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&COORD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END COORD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_adiabatic_dynamics_thermostat_slow_nose_force:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&FORCE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END FORCE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_adiabatic_dynamics_thermostat_slow_nose_mass:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&MASS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END MASS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_adiabatic_dynamics_thermostat_slow_nose_velocity:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&VELOCITY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END VELOCITY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_adiabatic_dynamics_thermostat_slow_nose:
    def __init__(self):
        self.params = {
                "LENGTH": None,
                "MTS": None,
                "TIMECON": None,
                "YOSHIDA": None,
                }
        self.status = False

        self.coord = cp2k_motion_md_adiabatic_dynamics_thermostat_slow_nose_coord()
        self.force = cp2k_motion_md_adiabatic_dynamics_thermostat_slow_nose_force()
        self.mass = cp2k_motion_md_adiabatic_dynamics_thermostat_slow_nose_mass()
        self.velocity = cp2k_motion_md_adiabatic_dynamics_thermostat_slow_nose_velocity()
        # basic setting
        self.params["LENGTH"] = 3
        self.params["MTS"] = 2
        self.params["TIMECON"] = 1.0e3
        self.params["YOSHIDA"] = 3

    def to_input(self, fout):
        fout.write("\t\t\t\t&NOSE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.coord.status == True:
            self.coord.to_input(fout)
        if self.force.status == True:
            self.force.to_input(fout)
        if self.mass.status == True:
            self.mass.to_input(fout)
        if self.velocity.status == True:
            self.velocity.to_input(fout)
        fout.write("\t\t\t\t&END NOSE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "VELOCITY":
                self.velocity.set_params({item: params[item]})
            elif item.split("-")[4] == "MASS":
                self.mass.set_params({item: params[item]})
            elif item.split("-")[4] == "FORCE":
                self.force.set_params({item: params[item]})
            elif item.split("-")[4] == "COORD":
                self.coord.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_adiabatic_dynamics_thermostat_slow:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.define_region = cp2k_motion_md_adiabatic_dynamics_thermostat_fast_define_region()
        self.nose = cp2k_motion_md_adiabatic_dynamics_thermostat_fast_nose()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&THERMOSTAT_SLOW\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.define_region.status == True:
            self.defin_region.to_input(fout)
        if self.nose.status == True:
            self.nose.to_input(fout)
        fout.write("\t\t\t&END THERMOSTAT_SLOW\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "DEFINE_REGION":
                self.define_region.set_params({item: params[item]})
            elif item.split("-")[3] == "NOSE":
                self.nose.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_adiabatic_dynamics:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.thermostat_fast = cp2k_motion_md_adiabatic_dynamics_thermostat_fast()
        self.thermostat_slow = cp2k_motion_md_adiabatic_dynamics_thermostat_slow()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t&ADIABATIC_DYNAMICS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.thermostat_fast.status == True:
            self.thermostat_fast.to_input(fout)
        if self.thermostat_slow.status == True:
            self.thermostat_slow.to_input(fout)
        fout.write("\t\t&END ADIABATIC_DYNAMICS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "THERMOSTAT_FAST":
                self.thermostat_fast.set_params({item: params[item]})
            elif item.split("-")[2] == "THERMOSTAT_SLOW":
                self.thermostat_slow.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_averages_print_averages_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
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


class cp2k_motion_md_averages_print_averages:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_averages_print_averages_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&PRINT_AVERAGES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END PRINT_AVERAGES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_md_averages_restart_averages:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&RESTART_AVERAGES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RESTART_AVERAGES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_averages:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.print_averages = cp2k_motion_md_averages_print_averages()
        self.restart_averages = cp2k_motion_md_averages_restart_averages()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t&AVERAGES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.print_averages.status == True:
            self.print_averages.to_input(fout)
        if self.restart_averages.status == True:
            self.restart_averages.to_input(fout)
        fout.write("\t\t&END AVERAGES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "PRINT_AVERAGES":
                self.print_averages.set_params({item: params[item]})
            elif item.split("-")[2] == "RESTART_AVERAGES":
                self.restart_averages.set_params({item: params[item]})
            else:
                pass




class cp2k_motion_md_barostat_mass:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&MASS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END MASS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_md_barostat_print_energy_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_barostat_print_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_barostat_print_energy_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_barostat_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.energy = cp2k_motion_md_barostat_print_energy()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.energy.status == True:
            self.energy.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "ENERGY":
                self.energy.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_barostat_thermostat_ad_langevin_chi:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&CHI\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END CHI\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_barostat_thermostat_ad_langevin_mass:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&MASS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END MASS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_barostat_thermostat_ad_langevin:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.chi = cp2k_motion_md_barostat_thermostat_ad_langevin_chi()
        self.mass = cp2k_motion_md_barostat_thermostat_ad_langevin_mass()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&AD_LANGEVIN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.chi.status == True:
            self.chi.to_input(fout)
        if self.mass.status == True:
            self.mass.to_input(fout)
        fout.write("\t\t\t\t&END AD_LANGEVIN\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CHI":
                self.chi.set_params({item: params[item]})
            elif item.split("-")[4] == "MASS":
                self.mass.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_barostat_thermostat_csvr_rng_init:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&RNG_INIT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END RNG_INIT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_barostat_thermostat_csvr_thermostat_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&THERMOSTAT_ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END THERMOSTAT_ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_barostat_thermostat_csvr:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.rng_init = cp2k_motion_md_barostat_thermostat_csvr_rng_init()
        self.thermostat_energy = cp2k_motion_md_barostat_thermostat_csvr_thermostat_energy()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&CSVR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.rng_init.status == True:
            self.rng_init.to_input(fout)
        if self.thermostat_energy.status == True:
            self.thermostat_energy.to_input(fout)
        fout.write("\t\t\t\t&END CSVR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "RNG_INIT":
                self.rng_init.set_params({item: params[item]})
            elif item.split("-")[4] == "THERMOSTAT_ENERGY":
                self.thermostat_energy.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_md_barostat_thermostat_gle_rng_init:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting
    
    def to_input(self, fout):
        fout.write("\t\t\t\t\t&RNG_INIT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END RNG_INIT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_barostat_thermostat_gle_s:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting
    
    def to_input(self, fout):
        fout.write("\t\t\t\t\t&S\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END S\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_barostat_thermostat_gle_thermostat_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting
    
    def to_input(self, fout):
        fout.write("\t\t\t\t\t&THERMOSTAT_ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END THERMOSTAT_ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_barostat_thermostat_gle:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.rng_init = cp2k_motion_md_barostat_thermostat_gle_rng_init()
        self.s = cp2k_motion_md_barostat_thermostat_gle_s()
        self.thermostat_energy = cp2k_motion_md_barostat_thermostat_gle_thermostat_energy()
        # basic setting
    
    def to_input(self, fout):
        fout.write("\t\t\t\t&GLE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.rng_init.status == True:
            self.rng_init.to_input(fout)
        if self.s.status == True:
            self.s.to_input(fout)
        if self.thermostat_energy.status == True:
            self.thermostat_energy.to_input(fout)
        fout.write("\t\t\t\t&END GLE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "RNG_INIT":
                self.rng_init.set_params({item: params[item]})
            elif item.split("-")[4] == "S":
                self.s.set_params({item: params[item]})
            elif item.split("-")[4] == "THERMOSTAT_ENERGY":
                self.thermostat_energy.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_md_barostat_thermostat_nose_coord:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&COORD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END COORD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_barostat_thermostat_nose_force:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&FORCE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END FORCE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_barostat_thermostat_nose_mass:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&MASS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END MASS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_barostat_thermostat_nose_velocity:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&VELOCITY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END VELOCITY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_barostat_thermostat_nose:
    def __init__(self):
        self.params = {
                "LENGTH": None,
                "MTS": None,
                "TIMECON": None,
                "YOSHIDA": None,
                }
        self.status = False

        self.coord = cp2k_motion_md_barostat_thermostat_nose_coord()
        self.force = cp2k_motion_md_barostat_thermostat_nose_force()
        self.mass = cp2k_motion_md_barostat_thermostat_nose_mass()
        self.velocity = cp2k_motion_md_barostat_thermostat_nose_velocity()
        # basic setting
        self.params["LENGTH"] = 3
        self.params["MTS"] = 2
        self.params["TIMECON"] = 1.0e3
        self.params["YOSHIDA"] = 3

    def to_input(self, fout):
        fout.write("\t\t\t\t&NOSE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.coord.status == True:
            self.coord.to_input(fout)
        if self.force.status == True:
            self.force.to_input(fout)
        if self.mass.status == True:
            self.mass.to_input(fout)
        if self.velocity.status == True:
            self.velocity.to_input(fout)
        fout.write("\t\t\t\t&END NOSE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "VELOCITY":
                self.velocity.set_params({item: params[item]})
            elif item.split("-")[4] == "MASS":
                self.mass.set_params({item: params[item]})
            elif item.split("-")[4] == "FORCE":
                self.force.set_params({item: params[item]})
            elif item.split("-")[4] == "COORD":
                self.coord.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_barostat_thermostat_print_energy_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_barostat_thermostat_print_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False
            
        self.each = cp2k_motion_md_barostat_thermostat_print_energy_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&ENERGY\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_md_barostat_thermostat_print_temperature_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_barostat_thermostat_print_temperature:
    def __init__(self):
        self.params = {
                }
        self.status = False
            
        self.each = cp2k_motion_md_barostat_thermostat_print_temperature_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&TEMPERATURE\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END TEMPERATURE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_barostat_thermostat_print_thermostat_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t&end EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_barostat_thermostat_print_thermostat_info:
    def __init__(self):
        self.params = {
                }
        self.status = False
            
        self.each = cp2k_motion_md_barostat_thermostat_print_thermostat_info_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&THERMOSTAT_INFO\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END THERMOSTAT_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_barostat_thermostat_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.energy = cp2k_motion_md_barostat_thermostat_print_energy()
        self.temperature = cp2k_motion_md_barostat_thermostat_print_temperature()
        self.thermostat_info = cp2k_motion_md_barostat_thermostat_print_thermostat_info()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.energy.status == True:
            self.energy.to_input(fout)
        if self.temperature.status == True:
            self.temperature.to_input(fout)
        if self.thermostat_info.status == True:
            self.thermostat_info.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "ENERGY":
                self.energy.set_params({item: params[item]})
            elif item.split("-")[4] == "TEMPERATURE":
                self.temperature.set_params({item: params[item]})
            elif item.split("-")[4] == "THERMOSTAT_INFO":
                self.thermostat_info.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_barostat_thermostat:
    def __init__(self):
        self.params = {
                "REGION": None,
                "TYPE": None,
                }
        self.status = False

        self.ad_langevin = cp2k_motion_md_barostat_thermostat_ad_langevin()
        self.csvr = cp2k_motion_md_barostat_thermostat_csvr()
        self.gle = cp2k_motion_md_barostat_thermostat_gle()
        self.nose = cp2k_motion_md_barostat_thermostat_nose()
        self.printout = cp2k_motion_md_barostat_thermostat_print()
        # basic setting
        self.params["TYPE"] = "NOSE"

    def to_input(self, fout):
        fout.write("\t\t\t&THERMOSTAT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.ad_langevin.status == True:
            self.add_langevin.to_input(fout)
        if self.csvr.status == True:
            self.csvr.to_input(fout)
        if self.gle.status == True:
            self.gle.to_input(fout)
        if self.nose.status == True:
            self.nose.to_input(fout)
        if self.params["TYPE"].upper() == "NOSE":
            self.nose.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)

        fout.write("\t\t\t&END THERMOSTAT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "ADD_LANGEVIN":
                self.add_langevin.set_params({item: params[item]})
            elif item.split("-")[3] == "CSVR":
                self.csvr.set_params({item: params[item]})
            elif item.split("-")[3] == "GLE":
                self.gle.set_params({item: params[item]})
            elif item.split("-")[3] == "NOSE":
                self.nose.set_params({item: params[item]})
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_md_barostat_velocity:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&VELOCITY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END VELOCITY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_barostat:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.mass = cp2k_motion_md_barostat_mass()
        self.printout = cp2k_motion_md_barostat_print()
        self.thermostat = cp2k_motion_md_barostat_thermostat()
        self.velocity = cp2k_motion_md_barostat_velocity()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t&BAROSTAT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.mass.status == True:
            self.mass.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.thermostat.status == True:
            self.thermostat.to_input(fout)
        if self.velocity.status == True:
            self.velocity.to_input(fout)
        fout.write("\t\t&END BAROSTAT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "MASS":
                self.mass.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[2] == "THERMOSTAT":
                self.thermostat.set_params({item: params[item]})
            elif item.split("-")[2] == "VELOCITY":
                self.velocity.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_cascade_atom_list:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&ATOM_LIST\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END ATOM_LIST\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_cascade:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.atom_list = cp2k_motion_md_cascade_atom_list()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t&CASCADE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.atom_list.status == True:
            self.atom_list.to_input(fout)
        fout.write("\t\t&END CASCADE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "ATOM_LIST":
                self.atom_list.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_md_initial_vibration:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t&INITIAL_VIBRATION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END INITIAL_VIBRATION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_langevin:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t&LANGEVIN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END LANGEVIN\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_msst:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t&MSST\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END MSST\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_print_center_of_mass_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_print_center_of_mass:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_print_center_of_mass_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&CENTER_OF_MASS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END CENTER_OF_MASS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_md_print_coefficients_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_print_coefficients:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_print_coefficients_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&COEFFICIENTS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END COEFFICIENTS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_md_print_energy_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_print_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_print_energy_each()
        # basic setting

    def to_input(self, fout):
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


class cp2k_motion_md_print_program_run_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_print_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_print_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_md_print_rotational_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_print_rotational_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_print_rotational_info_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&ROTATIONAL_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END ROTATIONAL_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_print_shell_energy_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_print_shell_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_print_shell_energy_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&SHELL_ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END SHELL_ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_md_print_temp_kind_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_print_temp_kind:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_print_temp_kind_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&TEMP_KIND\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END TEMP_KIND\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_print_temp_shell_kind_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_md_print_temp_shell_kind:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_print_temp_shell_kind_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&TEMP_SHELL_KIND\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END TEMP_SHELL_KIND\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.center_of_mass = cp2k_motion_md_print_center_of_mass()
        self.coefficients = cp2k_motion_md_print_coefficients()
        self.energy = cp2k_motion_md_print_energy()
        self.program_run_info = cp2k_motion_md_print_program_run_info()
        self.rotational_info = cp2k_motion_md_print_rotational_info()
        self.shell_energy = cp2k_motion_md_print_shell_energy()
        self.temp_kind = cp2k_motion_md_print_temp_kind()
        self.temp_shell_kind = cp2k_motion_md_print_temp_shell_kind()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.center_of_mass.status == True:
            self.center_of_mass.to_input(fout)
        if self.coefficients.status == True:
            self.coefficients.to_input(fout)
        if self.energy.status == True:
            self.energy.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.rotational_info.status == True:
            self.rotational_info.to_input(fout)
        if self.shell_energy.status == True:
            self.shell_energy.to_input(fout)
        if self.temp_kind.status == True:
            self.temp_kind.to_input(fout)
        if self.temp_shell_kind.status == True:
            self.temp_shell_kind.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "CENTER_OF_MASS":
                self.center_of_mass.set_params({item: params[item]})
            elif item.split("-")[2] == "COEFFICIENTS":
                self.coefficients.set_params({item: params[item]})
            elif item.split("-")[2] == "ENERGY":
                self.energy.set_params({item: params[item]})
            elif item.split("-")[2] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[2] == "ROTATIONAL_INFO":
                self.rotational_info.set_params({item: params[item]})
            elif item.split("-")[2] == "SHELL_ENENRGY":
                self.shell_energy.set_params({item: params[item]})
            elif item.split("-")[2] == "TEMP_KIND":
                self.temp_kind.set_params({item: params[item]})
            elif item.split("-")[2] == "TEMP_SHELL_KIND":
                self.temp_shell_kind.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_md_reftraj_msd_define_region:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t&DEFINE_REGION\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END DEFINE_REGION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_reftraj_msd:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.define_region = cp2k_motion_md_reftraj_msd_define_region()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&MSD\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.define_region.status == True:
            self.define_region.to_input(fout)
        fout.write("\t\t\t&END MSD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "DEFINE_REGION":
                self.define_region.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_md_reftraj_print_displaced_atom_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_reftraj_print_displaced_atom:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_reftraj_print_displaced_atom_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&DISPLACED_ATOM\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END DISPLACED_ATOM\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_reftraj_print_msd_kind_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_reftraj_print_msd_kind:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_reftraj_print_msd_kind_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&MSD_KIND\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END MSD_KIND\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_reftraj_print_msd_molecule_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_reftraj_print_msd_molecule:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_reftraj_print_msd_molecule_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&MSD_MOLECULE\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END MSD_MOLECULE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_md_reftraj_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.displaced_atom = cp2k_motion_md_reftraj_print_displaced_atom()
        self.msd_kind = cp2k_motion_md_reftraj_print_msd_kind()
        self.msd_molecule = cp2k_motion_md_reftraj_print_msd_molecule()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.displaced_atom.status == True:
            self.displaced_atom.to_input(fout)
        if self.msd_kind.status == True:
            self.msd_kind.to_input(fou)
        if self.msd_molecule.status == True:
            self.mds_molecul.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "DISPLACED_ATOM":
                self.displaced_atom.set_params({item: params[item]})
            elif item.split("_")[3] == "MSD_KIND":
                self.msd_kind.set_params({item: params[item]})
            elif item.split("-")[3] == "MSD_MOLECULE":
                self.msd_molecule.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_reftraj:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.msd = cp2k_motion_md_reftraj_msd()
        self.printout = cp2k_motion_md_reftraj_print()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t&reftraj\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.msd.status == True:
            self.msd.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t&end reftraj\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "MSD":
                self.msd.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_param({item: params[item]})
            else:
                pass


class cp2k_motion_md_respa:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t&RESPA\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END RESPA\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]






class cp2k_motion_md_shell_thermostat_ad_langevin_chi:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&CHI\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END CHI\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_shell_thermostat_ad_langevin_mass:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&MASS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END MASS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_shell_thermostat_ad_langevin:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.chi = cp2k_motion_md_shell_thermostat_ad_langevin_chi()
        self.mass = cp2k_motion_md_shell_thermostat_ad_langevin_mass()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&AD_LANGEVIN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.chi.status == True:
            self.chi.to_input(fout)
        if self.mass.status == True:
            self.mass.to_input(fout)
        fout.write("\t\t\t\t&END AD_LANGEVIN\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CHI":
                self.chi.set_params({item: params[item]})
            elif item.split("-")[4] == "MASS":
                self.mass.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_shell_thermostat_csvr_rng_init:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&RNG_INIT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END RNG_INIT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_shell_thermostat_csvr_thermostat_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&THERMOSTAT_ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END THERMOSTAT_ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_shell_thermostat_csvr:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.rng_init = cp2k_motion_md_shell_thermostat_csvr_rng_init()
        self.thermostat_energy = cp2k_motion_md_shell_thermostat_csvr_thermostat_energy()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&CSVR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.rng_init.status == True:
            self.rng_init.to_input(fout)
        if self.thermostat_energy.status == True:
            self.thermostat_energy.to_input(fout)
        fout.write("\t\t\t\t&END CSVR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "RNG_INIT":
                self.rng_init.set_params({item: params[item]})
            elif item.split("-")[4] == "THERMOSTAT_ENERGY":
                self.thermostat_energy.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_shell_thermostat_define_region:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&DEFINE_REGION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END DEFINE_REGION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_shell_thermostat_gle_rng_init:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting
    
    def to_input(self, fout):
        fout.write("\t\t\t\t\t&RNG_INIT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END RNG_INIT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_shell_thermostat_gle_s:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting
    
    def to_input(self, fout):
        fout.write("\t\t\t\t\t&S\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END S\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_shell_thermostat_gle_thermostat_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting
    
    def to_input(self, fout):
        fout.write("\t\t\t\t\t&THERMOSTAT_ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END THERMOSTAT_ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_shell_thermostat_gle:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.rng_init = cp2k_motion_md_shell_thermostat_gle_rng_init()
        self.s = cp2k_motion_md_shell_thermostat_gle_s()
        self.thermostat_energy = cp2k_motion_md_shell_thermostat_gle_thermostat_energy()
        # basic setting
    
    def to_input(self, fout):
        fout.write("\t\t\t\t&GLE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.rng_init.status == True:
            self.rng_init.to_input(fout)
        if self.s.status == True:
            self.s.to_input(fout)
        if self.thermostat_energy.status == True:
            self.thermostat_energy.to_input(fout)
        fout.write("\t\t\t\t&END GLE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "RNG_INIT":
                self.rng_init.set_params({item: params[item]})
            elif item.split("-")[4] == "S":
                self.s.set_params({item: params[item]})
            elif item.split("-")[4] == "THERMOSTAT_ENERGY":
                self.thermostat_energy.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_md_shell_thermostat_nose_coord:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&COORD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END COORD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_shell_thermostat_nose_force:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&FORCE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END FORCE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_shell_thermostat_nose_mass:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&MASS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END MASS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_shell_thermostat_nose_velocity:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&VELOCITY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END VELOCITY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_shell_thermostat_nose:
    def __init__(self):
        self.params = {
                "LENGTH": None,
                "MTS": None,
                "TIMECON": None,
                "YOSHIDA": None,
                }
        self.status = False

        self.coord = cp2k_motion_md_shell_thermostat_nose_coord()
        self.force = cp2k_motion_md_shell_thermostat_nose_force()
        self.mass = cp2k_motion_md_shell_thermostat_nose_mass()
        self.velocity = cp2k_motion_md_shell_thermostat_nose_velocity()
        # basic setting
        self.params["LENGTH"] = 3
        self.params["MTS"] = 2
        self.params["TIMECON"] = 1.0e3
        self.params["YOSHIDA"] = 3

    def to_input(self, fout):
        fout.write("\t\t\t\t&NOSE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.coord.status == True:
            self.coord.to_input(fout)
        if self.force.status == True:
            self.force.to_input(fout)
        if self.mass.status == True:
            self.mass.to_input(fout)
        if self.velocity.status == True:
            self.velocity.to_input(fout)
        fout.write("\t\t\t\t&END NOSE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "VELOCITY":
                self.velocity.set_params({item: params[item]})
            elif item.split("-")[4] == "MASS":
                self.mass.set_params({item: params[item]})
            elif item.split("-")[4] == "FORCE":
                self.force.set_params({item: params[item]})
            elif item.split("-")[4] == "COORD":
                self.coord.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_shell_thermostat_print_energy_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_shell_thermostat_print_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False
            
        self.each = cp2k_motion_md_shell_thermostat_print_energy_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&ENERGY\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_md_shell_thermostat_print_temperature_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_shell_thermostat_print_temperature:
    def __init__(self):
        self.params = {
                }
        self.status = False
            
        self.each = cp2k_motion_md_shell_thermostat_print_temperature_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&TEMPERATURE\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END TEMPERATURE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_shell_thermostat_print_thermostat_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_shell_thermostat_print_thermostat_info:
    def __init__(self):
        self.params = {
                }
        self.status = False
            
        self.each = cp2k_motion_md_shell_thermostat_print_thermostat_info_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&THERMOSTAT_INFO\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END THERMOSTAT_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_shell_thermostat_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.energy = cp2k_motion_md_shell_thermostat_print_energy()
        self.temperature = cp2k_motion_md_shell_thermostat_print_temperature()
        self.thermostat_info = cp2k_motion_md_shell_thermostat_print_thermostat_info()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.energy.status == True:
            self.energy.to_input(fout)
        if self.temperature.status == True:
            self.temperature.to_input(fout)
        if self.thermostat_info.status == True:
            self.thermostat_info.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "ENERGY":
                self.energy.set_params({item: params[item]})
            elif item.split("-")[4] == "TEMPERATURE":
                self.temperature.set_params({item: params[item]})
            elif item.split("-")[4] == "THERMOSTAT_INFO":
                self.thermostat_info.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_shell_thermostat:
    def __init__(self):
        self.params = {
                "REGION": None,
                "TYPE": None,
                }
        self.status = False

        self.ad_langevin = cp2k_motion_md_shell_thermostat_ad_langevin()
        self.csvr = cp2k_motion_md_shell_thermostat_csvr()
        self.define_region = cp2k_motion_md_shell_thermostat_define_region()
        self.gle = cp2k_motion_md_shell_thermostat_gle()
        self.nose = cp2k_motion_md_shell_thermostat_nose()
        self.printout = cp2k_motion_md_shell_thermostat_print()
        # basic setting
        self.params["TYPE"] = "NOSE"

    def to_input(self, fout):
        fout.write("\t\t&THERMOSTAT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.ad_langevin.status == True:
            self.add_langevin.to_input(fout)
        if self.csvr.status == True:
            self.csvr.to_input(fout)
        if self.define_region.status == True:
            self.define_region.to_input(fout)
        if self.gle.status == True:
            self.gle.to_input(fout)
        if self.nose.status == True:
            self.nose.to_input(fout)
        if self.params["TYPE"].upper() == "NOSE":
            self.nose.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)

        fout.write("\t\t&END THERMOSTAT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "ADD_LANGEVIN":
                self.add_langevin.set_params({item: params[item]})
            elif item.split("-")[3] == "CSVR":
                self.csvr.set_params({item: params[item]})
            elif item.split("-")[3] == "DEFINE_REGION":
                self.define_region.set_params({item: params[item]})
            elif item.split("-")[3] == "GLE":
                self.gle.set_params({item: params[item]})
            elif item.split("-")[3] == "NOSE":
                self.nose.set_params({item: params[item]})
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_shell:
    def __init__(self):
        self.params = {
                }
        self.status = False
            
        self.thermostat = cp2k_motion_md_shell_thermostat()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t&SHELL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.thermostat.status == True:
            self.thermostat.to_input(fout)
        fout.write("\t\t&END SHELL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "THERMOSTAT":
                self.thermostat.set_params({item: params[item]})


class cp2k_motion_md_thermal_region_define_region:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&DEFINE_REGION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END DEFINE_REGION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_thermal_region_print_langevin_regions_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
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

class cp2k_motion_md_thermal_region_print_langevin_regions:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_thermal_region_print_langevin_regions_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&LANGEVIN_REGIONS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END LANGEVIN_REGIONS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_md_thermal_region_print_temperature_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
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

class cp2k_motion_md_thermal_region_print_temperature:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_md_thermal_region_print_temperature_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&TEMPERATURE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END TEMPERATURE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_thermal_region_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.langevin_regions = cp2k_motion_md_thermal_region_print_langevin_regions()
        self.temperature = cp2k_motion_md_thermal_region_print_temperature()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.langevin_regions.status == True:
            self.langevin_region.to_input(fout)
        if self.temperature.status == True:
            self.temperature.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "LANGEVIN_REGIONS":
                self.langevin_regions.set_params({item: params[item]})
            elif item.split("-")[3] == "TEMPERATURE":
                self.temperature.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_thermal_region:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.define_region = cp2k_motion_md_thermal_region_define_region()
        self.printout = cp2k_motion_md_thermal_region_print()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t&THERMAL_REGION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.define_region.status == True:
            self.define_region.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t&END THERMAL_REGION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "DEFINE_REGION":
                self.define_region.set_params({item: params[item]})
            elif item.split("-")[2] == "PRITN":
                self.printout.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_md_thermostat_ad_langevin_chi:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&CHI\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END CHI\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_thermostat_ad_langevin_mass:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&MASS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END MASS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_thermostat_ad_langevin:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.chi = cp2k_motion_md_thermostat_ad_langevin_chi()
        self.mass = cp2k_motion_md_thermostat_ad_langevin_mass()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&AD_LANGEVIN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.chi.status == True:
            self.chi.to_input(fout)
        if self.mass.status == True:
            self.mass.to_input(fout)
        fout.write("\t\t\t&END AD_LANGEVIN\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "CHI":
                self.chi.set_params({item: params[item]})
            elif item.split("-")[3] == "MASS":
                self.mass.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_thermostat_csvr_rng_init:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&RNG_INIT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END RNG_INIT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_thermostat_csvr_thermostat_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&THERMOSTAT_ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END THERMOSTAT_ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_thermostat_csvr:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.rng_init = cp2k_motion_md_thermostat_csvr_rng_init()
        self.thermostat_energy = cp2k_motion_md_thermostat_csvr_thermostat_energy()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&CSVR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.rng_init.status == True:
            self.rng_init.to_input(fout)
        if self.thermostat_energy.status == True:
            self.thermostat_energy.to_input(fout)
        fout.write("\t\t\t&END CSVR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "RNG_INIT":
                self.rng_init.set_params({item: params[item]})
            elif item.split("-")[3] == "THERMOSTAT_ENERGY":
                self.thermostat_energy.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_thermostat_define_region:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&DEFINE_REGION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END DEFINE_REGION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_thermostat_gle_rng_init:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting
    
    def to_input(self, fout):
        fout.write("\t\t\t\t&RNG_INIT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END RNG_INIT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_thermostat_gle_s:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting
    
    def to_input(self, fout):
        fout.write("\t\t\t\t&S\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END S\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_thermostat_gle_thermostat_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting
    
    def to_input(self, fout):
        fout.write("\t\t\t\t&THERMOSTAT_ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END THERMOSTAT_ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_thermostat_gle:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.rng_init = cp2k_motion_md_thermostat_gle_rng_init()
        self.s = cp2k_motion_md_thermostat_gle_s()
        self.thermostat_energy = cp2k_motion_md_thermostat_gle_thermostat_energy()
        # basic setting
    
    def to_input(self, fout):
        fout.write("\t\t\t&GLE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.rng_init.status == True:
            self.rng_init.to_input(fout)
        if self.s.status == True:
            self.s.to_input(fout)
        if self.thermostat_energy.status == True:
            self.thermostat_energy.to_input(fout)
        fout.write("\t\t\t&END GLE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "RNG_INIT":
                self.rng_init.set_params({item: params[item]})
            elif item.split("-")[3] == "S":
                self.s.set_params({item: params[item]})
            elif item.split("-")[3] == "THERMOSTAT_ENERGY":
                self.thermostat_energy.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_md_thermostat_nose_coord:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&COORD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END COORD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_thermostat_nose_force:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&FORCE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END FORCE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_thermostat_nose_mass:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&MASS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END MASS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_thermostat_nose_velocity:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&VELOCITY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END VELOCITY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_md_thermostat_nose:
    def __init__(self):
        self.params = {
                "LENGTH": None,
                "MTS": None,
                "TIMECON": None,
                "YOSHIDA": None,
                }
        self.status = False

        self.coord = cp2k_motion_md_thermostat_nose_coord()
        self.force = cp2k_motion_md_thermostat_nose_force()
        self.mass = cp2k_motion_md_thermostat_nose_mass()
        self.velocity = cp2k_motion_md_thermostat_nose_velocity()
        # basic setting
        self.params["LENGTH"] = 3
        self.params["MTS"] = 2
        self.params["TIMECON"] = 1.0e3
        self.params["YOSHIDA"] = 3

    def to_input(self, fout):
        fout.write("\t\t\t&NOSE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.coord.status == True:
            self.coord.to_input(fout)
        if self.force.status == True:
            self.force.to_input(fout)
        if self.mass.status == True:
            self.mass.to_input(fout)
        if self.velocity.status == True:
            self.velocity.to_input(fout)
        fout.write("\t\t\t&END NOSE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "VELOCITY":
                self.velocity.set_params({item: params[item]})
            elif item.split("-")[3] == "MASS":
                self.mass.set_params({item: params[item]})
            elif item.split("-")[3] == "FORCE":
                self.force.set_params({item: params[item]})
            elif item.split("-")[3] == "COORD":
                self.coord.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_thermostat_print_energy_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_thermostat_print_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False
            
        self.each = cp2k_motion_md_thermostat_print_energy_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&energy\n")
        for item in self.params:
            if self.params[item] is not none:
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

class cp2k_motion_md_thermostat_print_temperature_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_thermostat_print_temperature:
    def __init__(self):
        self.params = {
                }
        self.status = False
            
        self.each = cp2k_motion_md_thermostat_print_temperature_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&TEMPERATURE\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END TEMPERATURE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_thermostat_print_thermostat_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&end EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_md_thermostat_print_thermostat_info:
    def __init__(self):
        self.params = {
                }
        self.status = False
            
        self.each = cp2k_motion_md_thermostat_print_thermostat_info_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&THERMOSTAT_INFO\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END THERMOSTAT_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_thermostat_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.energy = cp2k_motion_md_thermostat_print_energy()
        self.temperature = cp2k_motion_md_thermostat_print_temperature()
        self.thermostat_info = cp2k_motion_md_thermostat_print_thermostat_info()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.energy.status == True:
            self.energy.to_input(fout)
        if self.temperature.status == True:
            self.temperature.to_input(fout)
        if self.thermostat_info.status == True:
            self.thermostat_info.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "ENERGY":
                self.energy.set_params({item: params[item]})
            elif item.split("-")[3] == "TEMPERATURE":
                self.temperature.set_params({item: params[item]})
            elif item.split("-")[3] == "THERMOSTAT_INFO":
                self.thermostat_info.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_md_thermostat:
    def __init__(self):
        self.params = {
                "REGION": None,
                "TYPE": None,
                }
        self.status = False

        self.ad_langevin = cp2k_motion_md_thermostat_ad_langevin()
        self.csvr = cp2k_motion_md_thermostat_csvr()
        self.define_region = cp2k_motion_md_thermostat_define_region()
        self.gle = cp2k_motion_md_thermostat_gle()
        self.nose = cp2k_motion_md_thermostat_nose()
        self.printout = cp2k_motion_md_thermostat_print()
        # basic setting
        self.params["TYPE"] = "NOSE"

    def to_input(self, fout):
        fout.write("\t\t&THERMOSTAT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.ad_langevin.status == True:
            self.add_langevin.to_input(fout)
        if self.csvr.status == True:
            self.csvr.to_input(fout)
        if self.define_region.status == True:
            self.define_region.to_input(fout)
        if self.gle.status == True:
            self.gle.to_input(fout)
        if self.nose.status == True:
            self.nose.to_input(fout)
        if self.params["TYPE"].upper() == "NOSE":
            self.nose.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)

        fout.write("\t\t&END THERMOSTAT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "ADD_LANGEVIN":
                self.add_langevin.set_params({item: params[item]})
            elif item.split("-")[2] == "CSVR":
                self.csvr.set_params({item: params[item]})
            elif item.split("-")[2] == "DEFINE_REGION":
                self.define_region.set_params({item: params[item]})
            elif item.split("-")[2] == "GLE":
                self.gle.set_params({item: params[item]})
            elif item.split("-")[2] == "NOSE":
                self.nose.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_md_velocity_softening:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t&VELOCITY_SOFTENING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END VELOCITY_SOFTENING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]



class cp2k_motion_md:
    def __init__(self):
        self.params = {
                "ANGVEL_TOL": None,
                "ANGVEL_ZERO": None,
                "ANNEALING": None,
                "ANNEALING_CELL": None,
                "COMVEL_TOL": None,
                "DISPLACEMENT_TOL": None,
                "ECONS_START_VAL": None,
                "ENSEMBLE": None,
                "INITIAL_METHOD": None,
                "MAX_STEPS": None,
                "SCALE_TEMP_KIND": None,
                "STEPS": None,
                "STEP_START_VAL": None,
                "TEMPERATURE": None,
                "TEMPERATURE_ANNEALING": None,
                "TEMP_KIND": None,
                "TEMP_TOL": None,
                "TIMESTEP": None,
                "TIME_START_VAL": None,
                }
        self.status = False
        
        self.adiabatic_dynamics = cp2k_motion_md_adiabatic_dynamics()
        self.averages = cp2k_motion_md_averages()
        self.barostat = cp2k_motion_md_barostat()
        self.cascade = cp2k_motion_md_cascade()
        self.initial_vibration = cp2k_motion_md_initial_vibration()
        self.langevin = cp2k_motion_md_langevin()
        self.msst = cp2k_motion_md_msst()
        self.printout = cp2k_motion_md_print()
        self.reftraj = cp2k_motion_md_reftraj()
        self.respa = cp2k_motion_md_respa()
        self.shell = cp2k_motion_md_shell()
        self.thermal_region = cp2k_motion_md_thermal_region()
        self.thermostat = cp2k_motion_md_thermostat()
        self.velocity_softening = cp2k_motion_md_velocity_softening()
        # basic default setting
        self.params["TIMESTEP"] = 0.5
        self.params["STEPS"] = 1000

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&MD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, str(self.params[item])))
        if self.params["ENSEMBLE"] == "NVT":
            self.thermostat.to_input(fout)
        if self.params["ENSEMBLE"] == "REFTRAJ":
            self.reftraj.to_input(fout)
        if self.adiabatic_dynamics.status == True:
            self.adiabatic_dynamics.to_input(fout)
        if self.averages.status == True:
            self.averages.to_input(fout)
        if self.barostat.status == True:
            self.barostat.to_input(fout)
        if self.cascade.status == True:
            self.cascade.to_input(fout)
        if self.initial_vibration.status == True:
            self.initial_vibration.to_input(fout)
        if self.langevin.status == True:
            self.langevin.to_input(fout)
        if self.msst.status == True:
            self.msst.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.reftraj.status == True:
            self.reftraj.to_input(fout)
        if self.respa.status == True:
            self.respa.to_input(fout)
        if self.shell.status == True:
            self.shell.to_input(fout)
        if self.thermal_region.status == True:
            self.thermal_region.to_input(fout)
        if self.thermostat.status == True:
            self.thermostat.to_input(fout)
        if self.velocity_softening.status == True:
            self.velocity_softening.to_input(fout)
        fout.write("\t&END MD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "ADIABATIC_DYNAMICS":
                self.adiabatic_dynamics.set_params({item: params[item]})
            elif item.split("-")[1] == "AVERAGES":
                self.averages.set_params({item: params[item]})
            elif item.split("-")[1] == "BAROSTAT":
                self.barostat.set_params({item: params[item]})
            elif item.split("-")[1] == "CASCADE":
                self.cascade.set_params({item: params[item]})
            elif item.split("-")[1] == "INITIAL_VIBRATION":
                self.initial_vibration.set_params({item: params[item]})
            elif item.split("-")[1] == "LANGEVIN":
                self.langevin.set_params({item: params[item]})
            elif item.split("-")[1] == "MSST":
                self.msst.set_params({item: params[item]})
            elif item.split("-")[1] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[1] == "REFTRAJ":
                self.reftraj.set_params({item: params[item]})
            elif item.split("-")[1] == "RESPA":
                self.respa.set_params({item: params[item]})
            elif item.split("-")[1] == "SHELL":
                self.shell.set_params({item: params[item]})
            elif item.split("-")[1] == "THERMAL_REGION":
                self.thermal_region.set_params({item: params[item]})
            elif item.split("-")[1] == "THERMOSTAT":
                self.thermostat.set_params({item: params[item]})
            elif item.split("-")[1] == "VELOCITY_SOFTENING":
                self.velocity_softening.set_params({item: params[item]})
            else:
                pass
