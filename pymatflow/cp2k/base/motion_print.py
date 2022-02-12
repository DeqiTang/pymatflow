#!/usr/bin/env python
# _*_ coding: utf-8 _*_

class cp2k_motion_print_cell_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_print_cell:
    def __init__(self):
        self.params = {}
        self.status = False
        self.each = cp2k_motion_print_cell_each()
        
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&CELL\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&end CELL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_print_core_forces_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_print_core_forces:
    def __init__(self):
        self.params = {}
        self.status = False
        self.each = cp2k_motion_print_core_forces_each()
        
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&CORE_FORCES\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&end CORE_FORCES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_print_core_trajectory_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_print_core_trajectory:
    def __init__(self):
        self.params = {}
        self.status = False
        self.each = cp2k_motion_print_core_trajectory_each()
        
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&CORE_TRAJECTORY\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END CORE_TRAJECTORY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_print_core_velocities_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_print_core_velocities:
    def __init__(self):
        self.params = {}
        self.status = False
        self.each = cp2k_motion_print_core_velocities_each()
        
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&CORE_VELOCITIES\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END CORE_VELOCITIES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_print_forces_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_print_forces:
    def __init__(self):
        self.params = {}
        self.status = False
        self.each = cp2k_motion_print_forces_each()
        
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&FORCES\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END FORCES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_print_force_mixing_labels_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_print_force_mixing_labels:
    def __init__(self):
        self.params = {}
        self.status = False
        self.each = cp2k_motion_print_force_mixing_labels_each()
        
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&FORCE_MIXING_LABELS\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END FORCE_MIXING_LABELS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_print_mixed_energies_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_print_mixed_energies:
    def __init__(self):
        self.params = {}
        self.status = False
        self.each = cp2k_motion_print_mixed_energies_each()
        
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&MIXED_ENERGIES\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END MIXED_ENERGIES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_print_restart_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&each\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end each\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_print_restart:
    def __init__(self):
        self.params = {}
        self.status = False
        self.each = cp2k_motion_print_restart_each()
        
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&restart\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&end restart\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_print_restart_history_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&end EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_print_restart_history:
    def __init__(self):
        self.params = {}
        self.status = False

        self.each = cp2k_motion_print_restart_history_each()
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&restart_history\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&end restart_history\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})



class cp2k_motion_print_shell_forces_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_print_shell_forces:
    def __init__(self):
        self.params = {}
        self.status = False

        self.each = cp2k_motion_print_shell_forces_each()
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&SHELL_FORCES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END SHELL_FORCES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass




class cp2k_motion_print_shell_trajectory_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_print_shell_trajectory:
    def __init__(self):
        self.params = {}
        self.status = False

        self.each = cp2k_motion_print_shell_trajectory_each()
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&SHELL_TRAJECTORY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END SHELL_TRAJECTORY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_print_shell_velocities_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_print_shell_velocities:
    def __init__(self):
        self.params = {}
        self.status = False

        self.each = cp2k_motion_print_shell_velocities_each()
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&SHELL_VELOCITIES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END SHELL_VELOCITIES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_print_stress_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_print_stress:
    def __init__(self):
        self.params = {}
        self.status = False

        self.each = cp2k_motion_print_stress_each()
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&STRESS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END STRESS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_print_structure_data_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_print_structure_data:
    def __init__(self):
        self.params = {}
        self.status = False

        self.each = cp2k_motion_print_structure_data_each()
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&STRUCTURE_DATA\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END STRUCTURE_DATA\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_print_trajectory_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_print_trajectory:
    def __init__(self):
        self.params = {}
        self.status = False

        self.each = cp2k_motion_print_trajectory_each()
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&TRAJECTORY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END TRAJECTORY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass




class cp2k_motion_print_translation_vector_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_print_translation_vector:
    def __init__(self):
        self.params = {}
        self.status = False

        self.each = cp2k_motion_print_translation_vector_each()
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&TRANSLATION_VECTOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END TRANSLATION_VECTOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_print_velocities_each:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_print_velocities:
    def __init__(self):
        self.params = {}
        self.status = False

        self.each = cp2k_motion_print_velocities_each()
        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t&VELOCITIES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END VELOCITIES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_print:
    def __init__(self):
        self.status = False
        self.cell = cp2k_motion_print_cell()
        self.core_forces = cp2k_motion_print_core_forces()
        self.core_trajectory = cp2k_motion_print_core_trajectory()
        self.core_velocities = cp2k_motion_print_core_velocities()
        self.forces = cp2k_motion_print_forces()
        self.force_mixing_labels = cp2k_motion_print_force_mixing_labels()
        self.mixed_energies = cp2k_motion_print_mixed_energies()
        self.restart = cp2k_motion_print_restart()
        self.restart_history = cp2k_motion_print_restart_history()
        self.shell_forces = cp2k_motion_print_shell_forces()
        self.shell_trajectory = cp2k_motion_print_shell_trajectory()
        self.shell_velocities = cp2k_motion_print_shell_velocities()
        self.stress = cp2k_motion_print_stress()
        self.structure_data = cp2k_motion_print_structure_data()
        self.trajectory = cp2k_motion_print_trajectory()
        self.translation_vector = cp2k_motion_print_translation_vector()
        self.velocities = cp2k_motion_print_velocities()
        # basic setting
    
    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&PRINT\n")
        
        if self.cell.status == True:
            self.cell.to_input(fout)
        if self.core_forces.status == True:
            self.core_forces.to_input(fout)
        if self.core_trajectory.status == True:
            self.core_trajectory.to_input(fout)
        if self.core_velocities.status == True:
            self.core_velocities.to_input(fout)
        if self.forces.status == True:
            self.forces.to_input(fout)
        if self.force_mixing_labels.status == True:
            self.force_mixing_labels.to_input(fout)
        if self.mixed_energies.status == True:
            self.mixed_energies.to_input(fout)
        if self.restart.status == True:
            slef.restart.to_input(fout)
        if self.restart_history.status == True:
            self.restart_history.to_input(fout)
        if self.shell_forces.status == True:
            self.shell_forces.to_input(fout)
        if self.shell_trajectory.status == True:
            self.shell_trajectory.to_input(fout)
        if self.shell_velocities.status == True:
            self.shell_velocities.to_input(fout)
        if self.stress.status == True:
            self.stress.to_input(fout)
        if self.structure_data.status == True:
            self.structure_data.to_input(fout)
        if self.trajectory.status == True:
            self.trajectory.to_input(fout)
        if self.translation_vector.status == True:
            self.translation_vector.to_input(fout)
        if self.velocities.status == True:
            self.velocities.to_input(fout)
        fout.write("\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "CELL":
                self.cell.set_params({item: params[item]})
            elif item.split("-")[1] == "CORE_FORCES":
                self.core_forces.set_params({item: params[item]})
            elif item.split("-")[1] == "CORE_TRAJECTORY":
                self.core_trajectory.set_params({item: params[item]})
            elif item.split("-")[1] == "CORE_VELOCITIES":
                self.core_velocities.set_params({item: params[item]})
            elif item.split("-")[1] == "FORCES":
                self.forces.set_params({item: params[item]})
            elif item.split("-")[1] == "FORCE_MIXING_LABELS":
                self.force_mixing_labels.set_params({item: params[item]})
            elif item.split("-")[1] == "MIXED_ENERGIES":
                self.mixed_energies.set_params({item: params[item]})
            elif item.split("-")[1] == "RESTART":
                self.restart.set_params({item: params[item]})
            elif item.split("-")[1] == "RESTART_HISTORY":
                slef.restart_history.set_params({item: params[item]})
            elif item.split("-")[1] == "SHELL_FORCES":
                self.shell_forces.set_params({item: params[item]})
            elif item.split("-")[1] == "SHELL_TRAJECTORY":
                self.shell_trajectory.set_params({item: params[item]})
            elif item.split("-")[1] == "SHELL_VELOCITIES":
                slef.shell_velocities.set_params({item: params[item]})
            elif item.split("-")[1] == "STRESS":
                self.stress.set_params({item: params[item]})
            elif item.split("-")[1] == "STRUCTURE_DATA":
                self.structure_data.set_params({item: params[item]})
            elif item.split("-")[1] == "TRAJECTORY":
                self.trajectory.set_params({item: params[item]})
            elif item.split("-")[1] == "TRANSLATION_VECTOR":
                self.translation_vector.set_params({item: params[item]})
            elif item.split("-")[1] == "VELOCITIES":
                self.velocities.set_params({item: params[item]})
            else:
                pass
