#!/usr/bin/env python
# _*_ coding: utf-8 _*_



class cp2k_motion_constraint_collective_restraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RESTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RESTRAINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_constraint_collective:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.restraint = cp2k_motion_constraint_collective_restraint()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&COLLECTIVE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.restraint.status == True:
            self.restraint.to_input(fout)
        fout.write("\t\t&END COLLECTIVE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "RESTRAINT":
                self.restraint.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_constraint_colvar_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&COLVAR_RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END COLVAR_RESTART\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_constraint_constraint_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_constraint_constraint_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_constraint_constraint_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&CONSTRAINT_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END CONSTRAINT_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_constraint_fixed_atoms_restraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RESTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RESTRAINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_constraint_fixed_atoms:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.restraint = cp2k_motion_constraint_fixed_atoms_restraint()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&FIXED_ATOMS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.restraint.status == True:
            self.restraint.to_input(fout)
        fout.write("\t\t&END FIXED_ATOMS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "RESTRAINT":
                self.restraint.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_constraint_fix_atom_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&FIX_ATOM_RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END FIX_ATOM_RESTART\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_constraint_g3x3_restraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RESTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RESTRAINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_constraint_g3x3:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.restraint = cp2k_motion_constraint_g3x3_restraint()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&G3X3\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.restraint.status == True:
            self.restraint.to_input(fout)
        fout.write("\t\t&END G3X3\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "RESTRAINT":
                self.restraint.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_constraint_g4x6_restraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RESTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RESTRAINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_constraint_g4x6:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.restraint = cp2k_motion_constraint_g4x6_restraint()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&G4X6\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.restraint.status == True:
            self.restraint.to_input(fout)
        fout.write("\t\t&END G4X6\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "RESTRAINT":
                self.restraint.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_constraint_hbonds_restraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RESTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RESTRAINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_constraint_hbonds:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.restraint = cp2k_motion_constraint_hbonds_restraint()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&HBONDS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.restraint.status == True:
            self.restraint.to_input(fout)
        fout.write("\t\t&END HBONDS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "RESTRAINT":
                self.restraint.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_constraint_lagrange_multipliers_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_motion_constraint_lagrange_multipliers:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_constraint_lagrange_multipliers_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&LAGRANGE_MULTIPLIERS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END LAGRANGE_MULTIPLIERS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_constraint_virtual_site_restraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RESTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RESTRAINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_motion_constraint_virtual_site:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.restraint = cp2k_motion_constraint_virtual_site_restraint()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&VIRTUAL_SITE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.restraint.status == True:
            self.restraint.to_input(fout)
        fout.write("\t\t&END VIRTUAL_SITE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "RESTRAINT":
                self.restraint.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_constraint:
    def __init__(self):
        self.params = {
                "CONSTRAINT_INIT": None,
                "ROLL_TOLERANCE": None,
                "SHAKE_TOLERANCE": None,
                }
        self.status = False

        self.collective = cp2k_motion_constraint_collective()
        self.colvar_restart = cp2k_motion_constraint_colvar_restart()
        self.constraint_info = cp2k_motion_constraint_constraint_info()
        self.fixed_atoms = cp2k_motion_constraint_fixed_atoms()
        self.fix_atom_restart = cp2k_motion_constraint_fix_atom_restart()
        self.g3x3 = cp2k_motion_constraint_g3x3()
        self.g4x6 = cp2k_motion_constraint_g4x6()
        self.hbonds = cp2k_motion_constraint_hbonds()
        self.lagrange_multipliers = cp2k_motion_constraint_lagrange_multipliers()
        self.virtual_site = cp2k_motion_constraint_virtual_site()
            
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&CONSTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.collective.status == True:
            self.collective.to_input(fout)
        if self.colvar_restart.status == True:
            self.colvar_restart.to_input(fout)
        if self.constraint_info.status == True:
            self.constraitn_info.to_input(fout)
        if self.fixed_atoms.status == True:
            self.fixed_atoms.to_input(fout)
        if self.fix_atom_restart.status == True:
            self.fix_atom_restart.to_input(fout)
        if self.g3x3.status == True:
            self.g3x3.to_input(fout)
        if self.g4x6.status == True:
            self.g4x6.to_input(fout)
        if self.hbonds.status == True:
            self.hbonds.to_input(fout)
        if self.lagrange_multipliers.status == True:
            self.lagrange_multipliers.to_input(fout)
        if self.virtual_site.status == True:
            self.virtual_site.to_input(fout)
        fout.write("\t&END CONSTRAINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "COLLECTIVE":
                self.collective.set_params({item: params[item]})
            elif item.split("-")[1] == "COLVAR_RESTART":
                self.colvar_restart.set_params({item: params[item]})
            elif item.split("-")[1] == "CONSTRAINT_INFO":
                self.constraint_info.set_params({item: params[item]})
            elif item.split("-")[1] == "FIXED_ATOMS":
                self.fixed_atoms.set_params({item: params[item]})
            elif item.split("-")[1] == "FIX_ATOM_RESTART":
                self.fix_atom_restart.set_params({item: params[item]})
            elif item.split("-")[1] == "G3X3":
                self.g3x3.set_params({item: params[item]})
            elif item.split("-")[1] == "G4X6":
                self.g4x6.set_params({item: params[item]})
            elif item.split("-")[1] == "HBONDS":
                self.hbonds.set_params({item: params[item]})
            elif item.split("-")[1] == "LAGRANGE_MULTIPLIERS":
                self.lagrange_multipliers.set_params({item: params[item]})
            elif item.split("-")[1] == "VIRTUAL_SITE":
                self.virtual_site.set_params({item: params[item]})
            else:
                pass

