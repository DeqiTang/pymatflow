#!/usr/bin/env python
# _*_ coding: utf-8 _*_


class cp2k_properties_atomic:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        fout.write("\t\t&ATOMIC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END ATOMIC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_et_coupling_becke_constraint_a_atom_group:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        fout.write("\t\t\t\t&ATOM_GROUP\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END ATOM_GROUP\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_et_coupling_becke_constraint_a_dummy_atoms:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        fout.write("\t\t\t\t&DUMMY_ATOMS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END DUMMY_ATOMS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_et_coupling_becke_constraint_a_program_run_info_each:
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


class cp2k_properties_et_coupling_becke_constraint_a_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_et_coupling_becke_constraint_a_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
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


class cp2k_properties_et_coupling_becke_constraint_a:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.atom_group = cp2k_properties_et_coupling_becke_constraint_a_atom_group()
        self.dummy_atoms = cp2k_properties_et_coupling_becke_constraint_a_dummy_atoms()
        self.program_run_info = cp2k_properties_et_coupling_becke_constraint_a_program_run_info()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&BECKE_CONSTRAINT_A\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.atom_group.status == True:
            self.atom_group.to_input(fout)
        if self.dummy_atoms.status == True:
            self.dummy_atoms.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t\t&END BECKE_CONSTRAINT_A\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "ATOM_GROUP":
                self.atom_group.set_params({item: params[item]})
            elif item.split("-")[3] == "DUMMY_ATOMS":
                self.dummy_atoms.set_params({item: params[item]})
            elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_et_coupling_becke_constraint_b_atom_group:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        fout.write("\t\t\t\t&ATOM_GROUP\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END ATOM_GROUP\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_et_coupling_becke_constraint_b_dummy_atoms:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        fout.write("\t\t\t\t&DUMMY_ATOMS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END DUMMY_ATOMS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_et_coupling_becke_constraint_b_program_run_info_each:
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


class cp2k_properties_et_coupling_becke_constraint_b_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_et_coupling_becke_constraint_b_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
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


class cp2k_properties_et_coupling_becke_constraint_b:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.atom_group = cp2k_properties_et_coupling_becke_constraint_b_atom_group()
        self.dummy_atoms = cp2k_properties_et_coupling_becke_constraint_b_dummy_atoms()
        self.program_run_info = cp2k_properties_et_coupling_becke_constraint_b_program_run_info()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&BECKE_CONSTRAINT_B\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.atom_group.status == True:
            self.atom_group.to_input(fout)
        if self.dummy_atoms.status == True:
            self.dummy_atoms.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t\t&END BECKE_CONSTRAINT_B\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "ATOM_GROUP":
                self.atom_group.set_params({item: params[item]})
            elif item.split("-")[3] == "DUMMY_ATOMS":
                self.dummy_atoms.set_params({item: params[item]})
            elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_et_coupling_ddapc_restraint_a_program_run_info_each:
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


class cp2k_properties_et_coupling_ddapc_restraint_a_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_et_coupling_ddapc_restraint_a_program_run_info_each()
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
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


class cp2k_properties_et_coupling_ddapc_restraint_a:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.program_run_info = cp2k_properties_et_coupling_ddapc_restraint_a_program_run_info()
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t&DDAPC_RESTRAINT_A\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t\t&END DDAPC_RESTRAINT_A\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_et_coupling_ddapc_restraint_b_program_run_info_each:
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


class cp2k_properties_et_coupling_ddapc_restraint_b_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_et_coupling_ddapc_restraint_b_program_run_info_each()
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
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


class cp2k_properties_et_coupling_ddapc_restraint_b:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.program_run_info = cp2k_properties_et_coupling_ddapc_restraint_b_program_run_info()
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t&DDAPC_RESTRAINT_A\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t\t&END DDAPC_RESTRAINT_A\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_et_coupling_program_run_info_each:
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


class cp2k_properties_et_coupling_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_et_coupling_program_run_info_each()
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

class cp2k_properties_et_coupling_projection_block_print_mo_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_et_coupling_projection_block_print_mo_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_et_coupling_projection_block_print_mo_cubes_each()
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_et_coupling_projection_block_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.mo_cubes = cp2k_properties_et_coupling_projection_block_print_mo_cubes()
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.mo_cubes.status == True:
            self.mo_cubes.to_input(fout)
        fout.write("\t\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "MO_CUBES":
                slef.mo_cubes.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_et_coupling_projection_block:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.printout = cp2k_properties_et_coupling_projection_block_print()
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t&BLOCK\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t\t&END BLOCK\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_et_coupling_projection_program_run_info_each:
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


class cp2k_properties_et_coupling_projection_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_et_coupling_projection_program_run_info_each()
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
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

class cp2k_properties_et_coupling_projection:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.block = cp2k_properties_et_coupling_projection_block()
        self.program_run_info = cp2k_properties_et_coupling_projection_program_run_info()
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t&PROJECTION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.block.status == True:
            self.block.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t\t&END PROJECTION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "BLOCK":
                self.block.set_params({item: params[item]})
            elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_et_coupling:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.becke_constraint_a = cp2k_properties_et_coupling_becke_constraint_a()
        self.becke_constraint_b = cp2k_properties_et_coupling_becke_constraint_b()
        self.ddapc_restraint_a = cp2k_properties_et_coupling_ddapc_restraint_a()
        self.ddapc_restraint_b = cp2k_properties_et_coupling_ddapc_restraint_b()
        self.program_run_info = cp2k_properties_et_coupling_program_run_info()
        self.projection = cp2k_properties_et_coupling_projection()
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t&ET_COUPLING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.becke_constraint_a.status == True:
            slef.becke_constraint_a.to_input(fout)
        if self.becke_constraint_b.status == True:
            self.becke_constraint_b.to_input(fout)
        if self.ddapc_restraint_a.status == True:
            self.ddapc_restraint_a.to_input(fout)
        if self.ddapc_restraint_b.status == True:
            self.ddapc_restraint_b.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.projection.status == True:
            self.projection.to_input(fout)
        fout.write("\t\t&END ET_COUPLING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "BECKE_CONSTRAINT_A":
                self.becke_constraint_a.set_params({item: params[item]})
            elif item.split("-")[2] == "BECKE_CONSTRAINT_B":
                self.becke_constraint_b.set_params({item: params[item]})
            elif item.split("-")[2] == "DDAPC_RESTRAINT_A":
                self.ddapc_restraint_a.set_params({item: params[item]})
            elif item.split("-")[2] == "DDAPC_RESTRAINT_B":
                self.ddapc_restraint_b.set_params({item: params[item]})
            elif item.split("-")[2] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[2] == "PROJECTION":
                self.projection.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_fit_charge_each:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_fit_charge:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_fit_charge_each()

    def to_input(self, fout):
        fout.write("\t\t&FIT_CHARGE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END FIT_CHARGE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_resp_constraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&CONSTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END CONSTRAINT\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_resp_print_coord_fit_points_each:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        # basic setting


    def to_input(self, fout):
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


class cp2k_properties_resp_print_coord_fit_points:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_resp_print_coord_fit_points_each()
        
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t&COORD_FIT_POINTS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END COORD_FIT_POINTS\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "COORD_FIT_POINTS":
                self.coord_fit_points.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_resp_print_program_run_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        # basic setting


    def to_input(self, fout):
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


class cp2k_properties_resp_print_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_resp_print_program_run_info_each()
        
        # basic setting


    def to_input(self, fout):
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
            elif item.split("-")[4] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_resp_print_resp_charges_to_file_each:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        # basic setting


    def to_input(self, fout):
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

class cp2k_properties_resp_print_resp_charges_to_file:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_resp_print_resp_charges_to_file_each()
        
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t&RESP_CHARGES_TO_FILE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END RESP_CHARGES_TO_FILE\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_resp_print_v_resp_cube_each:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        # basic setting


    def to_input(self, fout):
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


class cp2k_properties_resp_print_v_resp_cube:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_resp_print_v_resp_cube_each()
        
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t&V_RESP_CUBE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END V_RESP_CUBE\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_resp_print:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.coord_fit_points = cp2k_properties_resp_print_coord_fit_points()
        self.program_run_info = cp2k_properties_resp_print_program_run_info()
        self.resp_charges_to_file = cp2k_properties_resp_print_resp_charges_to_file()
        self.v_resp_cube = cp2k_properties_resp_print_v_resp_cube()

        # basic setting
        self.v_resp_cube.status = True


    def to_input(self, fout):
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.coord_fit_points.status == True:
            self.coord_fit_points.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.resp_charges_to_file.status == True:
            self.resp_charges_to_file.to_input(fout)
        if self.v_resp_cube.status == True:
            self.v_resp_cube.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "COORD_FIT_POINTS":
                self.coord_fit_points.set_params({item: params[item]})
            elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[3] == "RESP_CHARGES_TO_FIEL":
                self.resp_charges_to_file.set_params({item: params[item]})
            elif item.split("-")[3] == "V_RESP_CUBE":
                self.v_resp_cube.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_resp_restraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&RESTARINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END RESTRAINT\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_resp_slab_sampling:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting
        self.params["RANGE"] = [0.3, 3.0]
        self.params["ATOM_LIST"] = [1]  # list like: [1, 2 ,3 , 4]
        self.params["SURF_DIRECTION"] = "Z"


    def to_input(self, fout):
        fout.write("\t\t\t&SLAB_SAMPLING\n")
        for item in self.params:
            if self.params[item] is not None:
                if item == "RANGE":
                    fout.write("\t\t\t\t%s %f %f\n" % (item, self.params[item][0], self.params[item][1]))
                elif item == "ATOM_LIST":
                    fout.write("\t\t\t\t%s " % item)
                    for atom in self.params[item]:
                        fout.write(" %d" % atom)
                    fout.write("\n")
                else:
                    fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END SLAB_SAMPLING\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_resp_sphere_sampling:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&SPHERE_SAMPLING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END SPHERE_SAMPLING\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_resp:
    def __init__(self):
        self.params = {
                "INTEGER_TOTAL_CHARGE": None,
                "RESTRAIN_HEAVIES_STRENGTH": None,
                "RESTRAIN_HEAVIES_TO_ZERO": None,
                "STRIDE": None,
                "USE_REPEAT_METHOD": None,
                "WIDTH": None,
                }
        self.status = False

        self.constraint = cp2k_properties_resp_constraint()
        self.printout = cp2k_properties_resp_print()
        self.restraint = cp2k_properties_resp_restraint()
        self.slab_sampling = cp2k_properties_resp_slab_sampling()
        self.sphere_sampling = cp2k_properties_resp_sphere_sampling()

        # basic setting
        self.printout.status = True
        self.slab_sampling.status = True

    def to_input(self, fout):
        fout.write("\t\t&RESP\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.constraint.status == True:
            self.constraint.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.restraint.status == True:
            self.restraint.to_input(fout)
        if self.slab_sampling.status == True:
            self.slab_sampling.to_input(fout)
        if self.sphere_sampling.status == True:
            self.sphere_samplint.to_input(fout)
        fout.write("\t\t&END RESP\n")

        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "CONSTRAINT":
                self.constraint.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[2] == "RESTRAINT":
                self.restraint.set_params({item: params[item]})
            elif item.split("-")[2] == "SLAB_SAMPLING":
                self.slab_sampling.set_params({item: params[item]})
            elif item.split("-")[2] == "SPHERE_SAMPLING":
                self.sphere_sampling.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_current_interpolator_conv_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_current_interpolator_conv_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_current_interpolator_conv_info_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&CONV_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END CONV_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_current_interpolator:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.conv_info = cp2k_properties_linres_current_interpolator_conv_info()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&INTERPOLATOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.conv_info.status == True:
            self.conv_info.to_input(fout)
        fout.write("\t\t\t\t&END INTERPOLATOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CONV_INFO":
                self.conv_info.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_current_print_current_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_current_print_current_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_current_print_current_cubes_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&CURRENT_CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END CURRENT_CUBES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_current_print_response_function_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_current_print_response_function_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_current_print_response_function_cubes_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&RESPONSE_FUNCTION_CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END RESPONSE_FUNCTION_CUBES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_current_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.current_cubes = cp2k_properties_linres_current_print_current_cubes()
        self.response_function_cubes = cp2k_properties_linres_current_print_response_function_cubes()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.current_cubes.status == True:
            self.current_cubes.to_input(fout)
        if self.response_function_cubes.status == True:
            self.response_function_cubes.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CURRENT_CUBES":
                self.current_cubes.set_params({item: params[item]})
            elif item.split("-")[4] == "RESPONSE_FUNCTION_CUBES":
                self.response_function_cubes.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_current:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.interpolator = cp2k_properties_linres_current_interpolator()
        self.printout = cp2k_properties_linres_current_print()
        
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&CURRENT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.interpolator.status == True:
            self.interpolator.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t&END CURRENT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "INTERPOLATOR":
                self.interpolator.set_params({item: params[item]})
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_interpolator_conv_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_epr_interpolator_conv_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_interpolator_conv_info_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&CONV_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END CONV_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_epr_interpolator:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.conv_info = cp2k_properties_linres_epr_interpolator_conv_info()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&INTERPOLATOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.conv_info.status == True:
            self.conv_info.to_input(fout)
        fout.write("\t\t\t\t&END INTERPOLATOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CONV_INFO":
                self.conv_info.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_adiabatic_rescaling:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&ADIABATIC_RESCALING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t\t&END ADIABATIC_RESCALING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_properties_linres_epr_print_g_tensor_xc_hf_hf_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 10:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_hf_hf_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_print_g_tensor_xc_hf_hf_info_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&HF_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t\t\t&END HF_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[8] == "EACH":
                self.each.set_params({item: params[item]})

class cp2k_properties_linres_epr_print_g_tensor_xc_hf_interaction_potential:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&INTERACTION_POTENTIAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t\t\t&END INTERACTION_POTENTIAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_hf_load_balance_print_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 11:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[10] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_hf_load_balance_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_print_g_tensor_xc_hf_load_balance_print_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 10:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[9] == "EACH":
                self.each.set_params({item: params[item]})

class cp2k_properties_linres_epr_print_g_tensor_xc_hf_load_balance:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.printout = cp2k_properties_linres_epr_print_g_tensor_xc_hf_load_balance_print()

        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&LOAD_BALANCE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t\t\t\t\t\t&END LOAD_BALANCE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[8] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_epr_print_g_tensor_xc_hf_memory:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&MEMORY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t\t\t&END MEMORY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_hf_periodic:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&PERIODIC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t\t\t&END PERIODIC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_hf_screening:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&SCREENING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t\t\t&END SCREENING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_properties_linres_epr_print_g_tensor_xc_hf:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.hf_info = cp2k_properties_linres_epr_print_g_tensor_xc_hf_hf_info()
        self.interaction_potential = cp2k_properties_linres_epr_print_g_tensor_xc_hf_interaction_potential()
        self.load_balance = cp2k_properties_linres_epr_print_g_tensor_xc_hf_load_balance()
        self.memory = cp2k_properties_linres_epr_print_g_tensor_xc_hf_memory()
        self.periodic = cp2k_properties_linres_epr_print_g_tensor_xc_hf_periodic()
        self.screening = cp2k_properties_linres_epr_print_g_tensor_xc_hf_screening()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&HF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.hf_info.status == True:
            self.hf_info.to_input(fout)
        if self.interaction_potential.status == True:
            self.interaction_potential.to_input(fout)
        if self.load_balance.status == True:
            self.load_balance.to_input(fout)
        if self.memory.status == True:
            self.memory.to_input(fout)
        if self.periodic.status == True:
            self.periodic.to_input(fout)
        if self.screening.status == True:
            self.screening.to_input(fout)
        fout.write("\t\t\t\t\t\t\t&END HF\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[7] == "HF_INFO":
                self.hf_info.set_params({item: parmas[item]})
            elif item.split("-")[7] == "INTERACTION_POTENTIAL":
                self.interaction_potential.set_params({item: params[item]})
            elif item.split("-")[7] == "LOAD_BALANCE":
                self.load_balance.set_params({item: params[item]})
            elif item.split("-")[7] == "MEMORY":
                self.memory.set_params({item: params[item]})
            elif item.split("-")[7] == "PERIODIC":
                self.periodic.set_params({item: params[item]})
            elif item.split("-")[7] == "SCREENING":
                self.screening.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_vdw_potential_pair_potential_print_dftd_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 11:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_vdw_potential_pair_potential_print_dftd:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_print_g_tensor_xc_vdw_potential_pair_potential_print_dftd_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t\t\t\t&PRINT_DFTD\n")
        for item in self.params:
            fout.write("\t\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t\t\t\t&END PRINT_DFTD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 10:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[9] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_vdw_potential_pair_potential:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.print_dftd = cp2k_properties_linres_epr_print_g_tensor_xc_vdw_potential_pair_potential_print_dftd()

        # basic setting
        self.params["PARAMETER_FILE_NAME"] = "dftd3.dat"

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t\t\t&PAIR_POTENTIAL\n")
        for item in self.params:
            fout.write("\t\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.print_dftd.status == True:
            self.print_dftd.to_input(fout)
        fout.write("\t\t\t\t\t\t\t\t&END PAIR_POTENTIAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[8] == "PRINT_DFTD":
                self.print_dftd.set_params({item: params[item]})

class cp2k_properties_linres_epr_print_g_tensor_xc_vdw_potential_non_local:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t\t\t&NON_LOCAL\n")
        for item in self.params:
            fout.write("\t\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t\t\t&END NON_LOCAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_properties_linres_epr_print_g_tensor_xc_vdw_potential:
    def __init__(self):
        self.params = {
                "POTENTIAL_TYPE": None,
                }
        self.status = False

        self.params["POTENTIAL_TYPE"] = "PAIR_POTENTIAL"
        self.pair_potential = cp2k_properties_linres_epr_print_g_tensor_xc_vdw_potential_pair_potential()
        self.non_local = cp2k_properties_linres_epr_print_g_tensor_xc_vdw_potential_non_local()

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t\t&VDW_POTENTIAL\n")
        for item in self.params:
            fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.params["POTENTIAL_TYPE"].upper() == "PAIR_POTENTIAL":
            self.pair_potential.to_input(fout)
        elif self.params["POTENTIAL_TYPE"].upper() == "NON_LOCAL":
            self.non_local.to_input(fout)
        fout.write("\t\t\t\t\t\t\t&END VDW_POTENTIAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[7] == "PAIR_POTENTIAL":
                self.pair_potential.set_params({item: params[item]})
            elif item.split("-")[7] == "NON_LOCAL":
                self.non_local.set_params({item: params[item]})



class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_cphf:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&CPHF\n")
        fout.write("\t\t\t\t\t\t\t\t&END CPHF\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_direct_canonical:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&WF_DIRECT_CANONICAL\n")
        fout.write("\t\t\t\t\t\t\t\t&END DIRECT_CANONICAL\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_eri_mme_cutoff_calib:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t&CUTOFF_CALIB\n")
        fout.write("\t\t\t\t\t\t\t\t\t&END CUTOFF_CALIB\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 10:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_eri_mme_eri_mme_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t\t&EACH\n")
        fout.write("\t\t\t\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 11:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_eri_mme_eri_mme_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_eri_mme_eri_mme_info_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t&ERI_MME_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t\t\t\t&END ERI_MME_INFO\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 10:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[9] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_eri_mme:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cutoff_calib = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_eri_mme_cutoff_calib()
        self.eri_mme_info = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_eri_mme_eri_mme_info()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&ERI_MME\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.cutoff_calib.status == True:
            self.cutoff_calib.to_input(fout)
        if self.eri_mme_info.status == True:
            self.eri_mme_info.to_input(fout)
        fout.write("\t\t\t\t\t\t\t\t&END ERI_MME\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[8] == "CUTOFF_CALIB":
                self.cutoff_calib.set_params({item: params[item]})
            elif item.split("-")[8] == "ERI_MME_INFO":
                self.eri_mme_info.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_interaction_potential:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&INTERACTION_POTENTIAL\n")
        fout.write("\t\t\t\t\t\t\t\t&END INTERACTION_POTENTIAL\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_mp2_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 10:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_mp2_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_mp2_info_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&MP2_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t\t\t&END MP2_INFO\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[8] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_opt_ri_basis:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&OPT_RI_BASIS\n")
        fout.write("\t\t\t\t\t\t\t\t&END OPT_RI_BASIS\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_laplace:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&RI_LAPLACE\n")
        fout.write("\t\t\t\t\t\t\t\t&END RI_LAPALACE\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_mp2:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&RI_MP2\n")
        fout.write("\t\t\t\t\t\t\t\t&END RI_MP2\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_hf_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t\t\t&EACH\n")
        fout.write("\t\t\t\t\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 12:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_hf_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_hf_info_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t\t&HF_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t\t\t\t\t&END HF_INFO\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 11:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[10] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_interaction_potential:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t\t&INTERACTION_POTENTIAL\n")
        fout.write("\t\t\t\t\t\t\t\t\t\t&END INTERACTION_POTENTIAL\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 11:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_load_balance_print_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t\t\t\t&EACH\n")
        fout.write("\t\t\t\t\t\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 13:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_load_balance_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_load_balance_print_each()

        # basic setting 
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t\t&LOAD_BALANCE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t\t\t\t\t&END LOAD_BALANCE\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 11:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[10] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_load_balance:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.printout = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_load_balance_print()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t\t&LOAD_BALANCE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t\t\t\t\t\t\t\t&END LOAD_BALANCE\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 11:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[10] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_memory:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t\t&MEMORY\n")
        fout.write("\t\t\t\t\t\t\t\t\t\t&END MEMORY\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 11:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_periodic:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t\t&PERIODIC\n")
        fout.write("\t\t\t\t\t\t\t\t\t\t&END PERIODIC\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 11:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_screening:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t\t&SCREENING\n")
        fout.write("\t\t\t\t\t\t\t\t\t\t&END SCREENING\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 11:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.hf_info = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_hf_info()
        self.interaction_potential = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_interaction_potential()
        self.load_balance = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_load_balance()
        self.memory = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_memory()
        self.periodic = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_periodic()
        self.screening = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf_screening()

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t&HF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.hf_info.status == True:    
            self.hf_info.to_input(fout)
        if self.interaction_potential.status == True:
            self.interaction_potential.to_input(fout)
        if self.load_balance.status == True:
            self.load_balance.to_input(fout)
        if self.memoery.status == True:
            self.memory.to_input(fout)
        if self.periodic.status == True:
            self.periodic.to_input(fout)
        if self.screening.status == True:
            self.screening.to_input()
        fout.write("\t\t\t\t\t\t\t\t\t&end HF\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 10:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[9] == "HF_INFO":
                self.hf_info.set_params({item: params[item]})
            elif item.split("-")[9] == "INTERACTION_POTENTIAL":
                self.interaction_potential.set_params({item: params[item]})
            elif item.split("-")[9] == "LOAD_BALANCE":
                self.load_balance.set_params({item: params[item]})
            elif item.split("-")[9] == "MEMORY":
                self.memory.set_params({item: params[item]})
            elif item.split("-")[9] == "PERIODIC":
                self.periodic.set_params({item: params[item]})
            elif item.split("-")[9] == "SCREENING":
                self.screening.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_im_time:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t&IM_TIME\n")
        fout.write("\t\t\t\t\t\t\t\t\t&end IM_TIME\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 10:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_ri_axk:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t&RI_AXK\n")
        fout.write("\t\t\t\t\t\t\t\t\t&end RI_AXK\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 10:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_ri_g0w0:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t\t&RI_G0W0\n")
        fout.write("\t\t\t\t\t\t\t\t\t&end RI_G0W0\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 10:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.hf = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_hf()
        self.im_time = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_im_time()
        self.ri_axk = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_ri_axk()
        self.ri_g0w0 = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa_ri_g0w0()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&RI_RPA\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.hf.status == True:
            self.hf.to_input(fout)
        if self.im_time.status == True:
            self.im_time.to_input(fout)
        if self.ri_axk.status == True:
            self.ri_axk.to_input(fout)
        if self.ri_g0w0.status == True:
            self.ri_g0w0.to_input(fout)
        fout.write("\t\t\t\t\t\t\t\t&end RI_RPA\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[8] == "HF":
                self.hf.set_params({item: params[item]})
            elif item.split("-")[8] == "IM_TIME":
                self.im_time.set_params({item: params[item]})
            elif item.split("-")[8] == "RI_AXK":
                self.ri_axk.set_params({item: params[item]})
            elif item.split("-")[8] == "RI_G0W0":
                self.ri_g0w0.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_wfc_gpw:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&WFC_GPW\n")
        fout.write("\t\t\t\t\t\t\t\t&END WFC_GPW\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cphf = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_cphf()
        self.direct_canonical = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_direct_canonical()
        self.eri_mme = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_eri_mme()
        self.interaction_potential = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_interaction_potential()
        self.mp2_info = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_mp2_info()
        self.opt_ri_basis = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_opt_ri_basis()
        self.ri_laplace = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_laplace()
        self.ri_mp2 = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_mp2()
        self.ri_rpa = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_ri_rpa()
        self.wfc_gpw = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation_wfc_gpw()

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&WF_CORRELATION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.cphf.status == True:
            self.cphf.to_input(fout)
        if self.direct_canonical.status == True:
            self.direct_canonical.to_input(fout)
        if self.eri_mme.status == True:
            self.eri_mme.to_input(fout)
        if self.interaction_potential.status == True:
            self.interaction_potential.to_input(fout)
        if self.mp2_info.status == True:
            self.mp2_info.to_input(fout)
        if self.opt_ri_basis.status == True:
            self.opt_ri_basis.to_input(fout)
        if self.ri_laplace.status == True:
            self.ri_laplace.to_iknput(fout)
        if self.ri_mp2.status == True:
            self.ri_mp2.to_input(fout)
        if self.ri_rpa.status == True:
            self.ri_rpa.to_input(fout)
        if self.wfc_gpw.status == True:
            self.wfc_gpw.to_input(fout)
        fout.write("\t\t\t\t\t\t\t&END WF_CORRELATION\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[7] == "CPHF":
                self.cphf.set_params({item: params[item]})
            elif item.split("-")[7] == "DIRECT_CANONICAL":
                self.direct_canonical.set_params({item: params[item]})
            elif item.split("-")[7] == "ERI_MME":
                self.eri_mme.set_params({item: params[item]})
            elif item.split("-")[7] == "INTERACTION_POTENTIAL":
                self.interaction_potential.set_params({item: params[item]})
            elif item.split("-")[7] == "MP2_INFO":
                self.mp2_info.set_params({item: params[item]})
            elif item.split("-")[7] == "OPT_RI_BASIS":
                self.opt_ri_basis.set_params({item: params[item]})
            elif item.split("-")[7] == "RI_LAPLACE":
                self.ri_laplace.set_params({item: params[item]})
            elif item.split("-")[7] == "RI_MP2":
                self.ri_mp2.set_params({item: params[item]})
            elif item.split("-")[7] == "RI_RPA":
                self.ri_rpa.set_params({item: params[item]})
            elif item.split("-")[7] == "WFC_GPW":
                self.wfc_gpw.set_params({item: params[item]})
            else:
                pass



class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_becke88:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&BECKE88\n")
        fout.write("\t\t\t\t\t\t\t\t&END BECKE88\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_becke88_lr:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&BECKE88_LR\n")
        fout.write("\t\t\t\t\t\t\t\t&END BECKE88_LR\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_becke88_lr_adiabatic:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&BECKE88_LR_ADIABATIC\n")
        fout.write("\t\t\t\t\t\t\t\t&END BECKE88_LR_ADIABATIC\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_becke97:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&BECKE97\n")
        fout.write("\t\t\t\t\t\t\t\t&END BECKE97\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_becke_roussel:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&BECKE_ROUSSEL\n")
        fout.write("\t\t\t\t\t\t\t\t&END BECKE_ROUSSEL\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_beef:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&BEEF\n")
        fout.write("\t\t\t\t\t\t\t\t&END BEEF\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_cs1:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&CS1\n")
        fout.write("\t\t\t\t\t\t\t\t&END CS1\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_gv09:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&GV09\n")
        fout.write("\t\t\t\t\t\t\t\t&END GV09\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_hcth:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&HCTH\n")
        fout.write("\t\t\t\t\t\t\t\t&END HCTH\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_ke_gga:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&KE_GGA\n")
        fout.write("\t\t\t\t\t\t\t\t&END KE_GGA\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_ke_libxc:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&KE_LIBXC\n")
        fout.write("\t\t\t\t\t\t\t\t&END KE_LIBXC\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_lda_hole_t_c_lr:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&LDA_HOLE_T_C_LR\n")
        fout.write("\t\t\t\t\t\t\t\t&END LDA_HOLE_T_C_LR\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_libxc:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&LIBXC\n")
        fout.write("\t\t\t\t\t\t\t\t&END LIBXC\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_lyp:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&LYP\n")
        fout.write("\t\t\t\t\t\t\t\t&END LYP\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_lyp_adiabatic:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&LYP_ADIABATIC\n")
        fout.write("\t\t\t\t\t\t\t\t&END LYP_ADIABATIC\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_optx:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&OPTX\n")
        fout.write("\t\t\t\t\t\t\t\t&END OPTX\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_p86c:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&P86C\n")
        fout.write("\t\t\t\t\t\t\t\t&END P86C\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_pade:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&PADE\n")
        fout.write("\t\t\t\t\t\t\t\t&END PADE\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_pbe:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&PBE\n")
        fout.write("\t\t\t\t\t\t\t\t&END PBE\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_pbe_hole_t_c_lr:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&PBE_HOLE_T_C_LR\n")
        fout.write("\t\t\t\t\t\t\t\t&END PBE_HOLE_T_C_LR\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_pw92:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&PW92\n")
        fout.write("\t\t\t\t\t\t\t\t&END PW92\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_pz81:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&PZ81\n")
        fout.write("\t\t\t\t\t\t\t\t&END PZ81\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_tf:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&TF\n")
        fout.write("\t\t\t\t\t\t\t\t&END TF\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_tfw:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&TFW\n")
        fout.write("\t\t\t\t\t\t\t\t&END TFW\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_tpss:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&TPSS\n")
        fout.write("\t\t\t\t\t\t\t\t&END TPSS\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_vwn:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&VWN\n")
        fout.write("\t\t\t\t\t\t\t\t&END VWN\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_xalpha:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&XALPHA\n")
        fout.write("\t\t\t\t\t\t\t\t&END XALPHA\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_xgga:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&XGGA\n")
        fout.write("\t\t\t\t\t\t\t\t&END XGGA\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_xwpbe:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&XWPBE\n")
        fout.write("\t\t\t\t\t\t\t\t&END XWPBE\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional:
    def __init__(self):
        self.section = "PBE"
        self.params = {
                }
        self.status = False

        self.becke88 =  cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_becke88()
        self.becke88_lr = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_becke88_lr()
        self.becke88_lr_adiabatic = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_becke88_lr_adiabatic()
        self.becke97 = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_becke97()
        self.becke_roussel = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_becke_roussel()
        self.beef = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_beef()
        self.cs1 = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_cs1()
        self.gv09 = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_gv09()
        self.hcth = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_hcth()
        self.ke_gga = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_ke_gga()
        self.ke_libxc = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_ke_libxc()
        self.lda_hole_t_c_lr = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_lda_hole_t_c_lr()
        self.libxc = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_libxc()
        self.lyp = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_lyp()
        self.lyp_adiabatic = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_lyp_adiabatic()
        self.optx = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_optx()
        self.p86c = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_p86c()
        self.pade = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_pade()
        self.pbe = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_pbe()
        self.pbe_hole_t_c_lr = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_pbe_hole_t_c_lr()
        self.pw92 = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_pw92()
        self.pz81 = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_pz81()
        self.tf = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_tf()
        self.tfw = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_tfw()
        self.tpss = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_tpss()
        self.vwn = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_vwn()
        self.xalpha = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_xalpha()
        self.xgga = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_xgga()
        self.xwpbe = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional_xwpbe()

        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&XC_FUNCTIONAL %s\n" % self.section)
        
        if self.becke88.status == True:
            self.becke88.to_input(fout)
        if self.becke88_lr == True:
            self.becke88_lr.to_input(fout)
        if self.becke88_lr_adiabatic == True:
            self.becke88_lr_adiabatic.to_input(fout)
        if self.becke97 == True:
            self.becke97.to_input(fout)
        if self.becke_roussel == True:
            self.becke_roussel.to_input(fout)
        if self.beef == True:
            self.beef.to_input(fout)
        if self.cs1 == True:
            self.cs1.to_input(fout)
        if self.gv09 == True:
            self.gv09.to_input(fout)
        if self.hcth == True:
            self.hcth.to_input(fout)
        if self.ke_gga == True:
            self.ke_gga.to_input(fout)
        if self.ke_libxc == True:
            self.ke_libxc.to_input(fout)
        if self.lda_hole_t_c_lr == True:
            self.lda_hole_t_c_lr.to_input(fout)
        if self.libxc == True:
            self.libxc.to_input(fout)
        if self.lyp == True:
            self.lyp.to_input(fout)
        if self.lyp_adiabatic == True:
            self.lyp_adiabatic.to_input(fout)
        if self.optx == True:
            self.optx.to_input(fout)
        if self.p86c == True:
            self.p86c.to_input(fout)
        if self.pade == True:
            self.pade.to_input(fout)
        if self.pbe == True:
            self.pbe.to_input(fout)
        if self.pbe_hole_t_c_lr == True:
            self.pbe_hole_t_c_lr.to_input(fout)
        if self.pw92 ==True:
            self.pw92.to_input(fout)
        if self.pz81 == True:
            self.pz81.to_input(fout)
        if self.tf == True:
            self.tf.to_input(fout)
        if self.tfw == True:
            self.tfw.to_input(fout)
        if self.tpss == True:
            self.tpss.to_input(fout)
        if self.vwn == True:
            self.vwn.to_input(fout)
        if self.xalpha == True:
            self.xalpha.to_input(fout)
        if self.xgga == True:
            self.xgga.to_input(fout)
        if self.xwpbe == True:
            self.xwpbe.to_input(fout)

        fout.write("\t\t\t\t\t\t\t&END XC_FUNCTIONAL\n")

    def set_params(self, params):
        """
        set_params for xc_functional is different from many other
        set_params, as it deal with the key 'XC_FUNCTIONAL' only
        """
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[7] == "BECKE88":
                self.becke88.set_params({item: params[item]})
            elif item.split("-")[7] == "BECKE88_LR":
                self.becke88_lr.set_params({item: params[item]})
            elif item.split("-")[7] == "BECKE88_LR_ADIABATIC":
                self.becke88_lr_adiabatic.set_params({item: params[item]})
            elif item.split("-")[7] == "BECKE97":
                self.becke97.set_params({item: params[item]})
            elif item.split("-")[7] == "BECKE_ROUSSEL":
                self.becke_roussel.set_params({item: params[item]})
            elif item.split("-")[7] == "BEEF":
                self.beef.set_params({item: params[item]})
            elif item.split("-")[7] == "CS1":
                self.cs1.set_params({item: params[item]})
            elif item.split("-")[7] == "GV09":
                self.gv09.set_params({item: params[item]})
            elif item.split("-")[7] == "HCTH":
                self.hcth.set_params({item: params[item]})
            elif item.split("-")[7] == "KE_GGA":
                self.ke_gga.set_params({item: params[item]})
            elif item.split("-")[7] == "KE_LIBXC":
                self.ke_libxc.set_params({item: params[item]})
            elif item.split("-")[7] == "LDA_HOLE_T_C_LR":
                self.lda_hole_t_c_lr.set_params({item: params[item]})
            elif item.split("-")[7] == "LIBXC":
                self.libxc.set_params({item: params[item]})
            elif item.split("-")[7] == "LYP":
                self.lyp.set_params({item: params[item]})
            elif item.split("-")[7] == "LYP_ADIABATIC":
                self.lyp_adiabatic.set_params({item: params[item]})
            elif item.split("-")[7] == "OPTX":
                self.optx.set_params({item: params[item]})
            elif item.split("-")[7] == "P86C":
                self.p86c.set_params({item: params[item]})
            elif item.split("-")[7] == "PADE":
                self.pade.set_params({item: params[item]})
            elif item.split("-")[7] == "PBE":
                self.pbe.set_params({item: params[item]})
            elif item.split("-")[7] == "PBE_HOLE_T_C_LR":
                self.pbe_hole_t_c_lr.set_params({item: params[item]})
            elif item.split("-")[7] == "PW92":
                self.pw92.set_params({item: params[item]})
            elif item.split("-")[7] == "PZ81":
                self.pz81.set_params({item: params[item]})
            elif item.split("-")[7] == "TF":
                self.tf.set_params({item: params[item]})
            elif item.split("-")[7] == "TFW":
                self.tfw.set_params({item: params[item]})
            elif item.split("-")[7] == "TPSS":
                self.tpss.set_params({item: params[item]})
            elif item.split("-")[7] == "VWN":
                self.vwn.set_params({item: params[item]})
            elif item.split("-")[7] == "XALPHA":
                self.xalpha.set_params({item: params[item]})
            elif item.split("-")[7] == "XGGA":
                self.xgga.set_params({item: params[item]})
            elif item.split("-")[7] == "XWPBE":
                self.xwpbe.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_grid:
    def __init__(self):
        self.params = {
                "XC_DERIV": "NN10_SMOOTH",
                "XC_SMOOTH_RHO": "NN10",
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&XC_GRID\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t\t&END XC_GRID\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_properties_linres_epr_print_g_tensor_xc_xc_potential_saop:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&SAOP\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t\t\t&END SAOP\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_properties_linres_epr_print_g_tensor_xc_xc_potential:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.saop = cp2k_properties_linres_epr_print_g_tensor_xc_xc_potential_saop()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&XC_POTENTIAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.saop.status == True:
            self.saop.to_input(fout)
        fout.write("\t\t\t\t\t\t\t&END XC_POTENTIAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[7] == "SAOP":
                self.saop.set_params({item: params[item]})


class cp2k_properties_linres_epr_print_g_tensor_xc:
    def __init__(self):
        self.params = {
                "DENSITY_CUTOFF": None,
                "DENSITY_SMOOTH_CUTOFF_RANGE": None,
                "FUNCTIONAL_ROUTINE": None,
                "GRADIENT_CUTOFF": None,
                "TAU_CUTOFF": None,
                }
        self.status = False

        self.adiabatic_rescaling = cp2k_properties_linres_epr_print_g_tensor_xc_adiabatic_rescaling()
        self.hf = cp2k_properties_linres_epr_print_g_tensor_xc_hf()
        self.vdw_potential = cp2k_properties_linres_epr_print_g_tensor_xc_vdw_potential()
        self.wf_correlation = cp2k_properties_linres_epr_print_g_tensor_xc_wf_correlation()
        self.xc_functional = cp2k_properties_linres_epr_print_g_tensor_xc_xc_functional()
        self.xc_grid = cp2k_properties_linres_epr_print_g_tensor_xc_xc_grid()
        self.xc_potential = cp2k_properties_linres_epr_print_g_tensor_xc_xc_potential()

        # basic setting
        self.xc_functional.status = True
        self.vdw_potential.status = False
        self.xc_grid.status = False
        self.adiabatic_rescaling.status = False
        self.hf.status = False
        self.wf_correlation.status = False
        self.xc_potential.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&XC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s" % (item, str(self.params[item])))
        if self.adiabatic_rescaling.status == True:
            self.adiabatic_rescaling.to_input(fout)
        if self.hf.status == True:
            self.hf.to_input(fout)
        if self.xc_functional.status == True:
            self.xc_functional.to_input(fout)
        if self.vdw_potential.status == True:
            self.vdw_potential.to_input(fout)
        if self.wf_correlation.status == True:
            self.wf_correlation.to_input(fout)
        if self.xc_grid.status == True:
            self.xc_grid.to_input(fout)
        if self.xc_potential.status == True:
            self.xc_potential.to_input(fout)
        fout.write("\t\t\t\t\t\t&END XC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                if item.split("-")[-1] == "XC_FUNCTIONAL":
                    self.xc_functional.section = params[item]
                elif item.split("-")[-1] == "VDW_POTENTIAL":
                    self.vdw_potential.section = params[item]
                else:
                    self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[6] == "XC_FUNCTIONAL" and len(item.split("-")) > 3:
                self.xc_functional.set_params({item: params[item]})
            elif item.split("-")[6] == "VDW_POTENTIAL":
                self.vdw_potential.set_params({item: params[item]})
            elif item.split("-")[6] == "ADIABATIC_RESCALING":
                self.adiabatic_rescaling.set_params({item: params[item]})
            elif item.split("-")[6] == "HF":
                self.hf.set_params({item: params[item]})
            elif item.split("-")[6] == "WF_CORRELATION":
                self.wf_correlation.set_params({item: params[item]})
            elif item.split("-")[6] == "XC_GRID":
                self.xc_grid.set_params({item: params[item]})
            elif item.split("-")[6] == "XC_POTENTIAL":
                self.xc_potential.set_params({item: params[item]})
            else:
                pass



class cp2k_properties_linres_epr_print_g_tensor_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_epr_print_g_tensor:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_print_g_tensor_each()
        self.xc = cp2k_properties_linres_epr_print_g_tensor_xc()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&G_TENSOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        if self.xc.status == True:
            self.xc.to_input(fout)
        fout.write("\t\t\t\t\t&END G_TENSOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            elif item.split("-")[5] == "XC":
                self.xc.set_params({item: params[item]})
            else:
                pass



class cp2k_properties_linres_epr_print_nablavks_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_epr_print_nablavks_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_print_nablavks_cubes_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&NABLAVKS_CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END NABLAVKS_CUBES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_print_response_function_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_epr_print_response_function_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_print_response_function_cubes_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&RESPONSE_FUNCTION_CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END RESPONSE_FUNCTION_CUBES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_epr_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.g_tensor = cp2k_properties_linres_epr_print_g_tensor()
        self.nablavks_cubes = cp2k_properties_linres_epr_print_nablavks_cubes()
        self.response_function_cubes = cp2k_properties_linres_epr_print_response_function_cubes()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.g_tensor.status == True:
            self.g_tensor.to_input(fout)
        if self.nablavks_cubes.status == True:
            self.nablavks_cubes.to_input(fout)
        if self.response_function_cubes.status == True:
            self.response_function_cubes.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "G_TENSOR":
                self.g_tensor.set_params({item: params[item]})
            elif item.split("-")[4] == "NABLAVKS_CUBES":
                self.nablavks_cubes.set_params({item: params[item]})
            elif item.split("-")[4] == "RESPONSE_FUNCTION_CUBES":
                self.response_function_cubes.set_params({item: params[item]})
            else:
                pass



class cp2k_properties_linres_epr:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.interpolator  = cp2k_properties_linres_epr_interpolator()
        self.printout = cp2k_properties_linres_epr_print()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&EPR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.interpolator.status == True:
            self.interpolator.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t&END EPR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "INTERPOLATOR":
                self.interpolator.set_params({item: params[item]})
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_localize_print_loc_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&LOC_RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END LOC_RESTART\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_localize_print_loc_restart_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_localize_print_loc_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_localize_print_loc_restart_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&LOC_RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END LOC_RESTART\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_localize_print_molecular_dipoles_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_localize_print_molecular_dipoles:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_localize_print_molecular_dipoles_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&MOLECULAR_DIPOLES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END MOLECULAR_DIPOLES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_localize_print_molecular_states_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_localize_print_molecular_states_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_localize_print_molecular_states_cubes_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t&END CUBES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[6] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_localize_print_molecular_states_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_localize_print_molecular_states:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_localize_print_molecular_states_each()
        self.cubes = cp2k_properties_linres_localize_print_molecular_states_cubes()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&MOLECULAR_STATES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        if self.cubes.status == True:
            self.cubes.to_input(fout)
        fout.write("\t\t\t\t\t&END MOLECULAR_STATES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            elif item.split("-")[5] == "CUBES":
                self.cubes.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_localize_print_program_run_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_localize_print_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_localize_print_program_run_info_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_localize_print_total_dipole_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_localize_print_total_dipole:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_localize_print_total_dipole_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&TOTAL_DIPOLE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END TOTAL_DIPOLE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_localize_print_wannier_centers_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_localize_print_wannier_centers:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_localize_print_wannier_centers_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&WANNIER_CENTERS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END WANNIER_CENTERS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_localize_print_wannier_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_localize_print_wannier_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_localize_print_wannier_cubes_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&WANNIER_CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END WANNIER_CUBES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_localize_print_wannier_spreads_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_localize_print_wannier_spreads:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_localize_print_wannier_spreads_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&WANNIER_SPREADS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END WANNIER_SPREADS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_localize_print_wannier_states_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_localize_print_wannier_states_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_localize_print_wannier_states_cubes_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t&END CUBES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[6] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_localize_print_wannier_states_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_localize_print_wannier_states:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_localize_print_wannier_states_each()
        self.cubes = cp2k_properties_linres_localize_print_wannier_states_cubes()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&WANNIER_STATES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.cubes.status == True:
            self.cubes.to_input(fout)
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END WANNIER_STATES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "CUBES":
                self.cubes.set_params({item: params[item]})
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_localize_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.loc_restart = cp2k_properties_linres_localize_print_loc_restart()
        self.molecular_states = cp2k_properties_linres_localize_print_molecular_states()
        self.program_run_info = cp2k_properties_linres_localize_print_program_run_info()
        self.total_dipole = cp2k_properties_linres_localize_print_total_dipole()
        self.wannier_centers = cp2k_properties_linres_localize_print_wannier_centers()
        self.wannier_cubes = cp2k_properties_linres_localize_print_wannier_cubes()
        self.wannier_spreads = cp2k_properties_linres_localize_print_wannier_spreads()
        self.wannier_states = cp2k_properties_linres_localize_print_wannier_states()

        # basic setting

    def to_input(self, fout):
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
            elif item.split("-")[4] == "WANNIER_CUEBS":
                self.wannier_cubes.set_params({item: params[item]})
            elif item.split("-")[4] == "WANNIER_SPREADS":
                self.wannier_spreads.set_params({item: params[item]})
            elif item.split("-")[4] == "WANNIER_STATES":
                self.wannier_states.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_localize:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.printout = cp2k_properties_linres_localize_print()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&LOCALIZE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t&END LOCALIZE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_nmr_interpolator_conv_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_nmr_interpolator_conv_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_nmr_interpolator_conv_info_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&CONV_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END CONV_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_nmr_interpolator:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.conv_info = cp2k_properties_linres_nmr_interpolator_conv_info()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&INTERPOLATOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.conv_info.status == True:
            self.conv_info.to_input(fout)
        fout.write("\t\t\t\t&END INTERPOLATOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CONV_INFO":
                self.conv_info.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_nmr_print_chi_tensor_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_nmr_print_chi_tensor:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_nmr_print_chi_tensor_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&CHI_TENSOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END CHI_TENSOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_nmr_print_response_function_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_nmr_print_response_function_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_nmr_print_response_function_cubes_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&RESPONSE_FUNCTION_CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END RESPONSE_FUNCTION_CUBES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_nmr_print_shielding_tensor_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_nmr_print_shielding_tensor:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_nmr_print_shielding_tensor_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&SHIELDING_TENSOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END SHIELDING_TENSOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_nmr_print:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.chi_tensor = cp2k_properties_linres_nmr_print_response_function_cubes()
        self.response_function_cubes = cp2k_properties_linres_nmr_print_response_function_cubes()
        self.shielding_tensor = cp2k_properties_linres_nmr_print_shielding_tensor()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.chi_tensor.status == True:
            self.chi_tensor.to_input(fout)
        if self.response_function_cubes.status == True:
            self.response_function_cubes.to_input(fout)
        if self.shielding_ternsor.status == True:
            self.shielding_tensor.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CHI_TENSOR":
                self.chi_tensor.set_params({item: params[item]})
            elif item.split("-")[4] == "RESPONSE_FUNCTION_CUBES":
                self.response_function_cubes.set_params({item: params[item]})
            elif item.split("-")[4] == "SHIELDING_TENSOR":
                self.shielding_tensor.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_nmr:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.interpolator = cp2k_properties_linres_nmr_interpolator
        self.printout =cp2k_properties_linres_nmr_print()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&NMR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.interpolator.status == True:
            self.interpolator.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t&END NMR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "INTERPOLATOR":
                self.interpolator.set_params({item: params[item]})
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_polar_interpolator_conv_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_polar_interpolator_conv_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_polar_interpolator_conv_info_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&CONV_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END CONV_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_polar_interpolator:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.conv_info = cp2k_properties_linres_polar_interpolator_conv_info()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&INTERPOLATOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.conv_info.status == True:
            self.conv_info.to_input(fout)
        fout.write("\t\t\t\t&END INTERPOLATOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CONV_INFO":
                self.conv_info.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_polar_print_polar_matrix_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_polar_print_polar_matrix:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_polar_print_polar_matrix_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&POLAR_MATRIX\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END POLAR_MATRIX\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_polar_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.polar_matrix = cp2k_properties_linres_polar_print_polar_matrix()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.polar_matrix.status == True:
            self.polar_matrix.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "POLAR_MATRIX":
                self.polar_matrix.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_polar:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.interpolator = cp2k_properties_linres_polar_interpolator()
        self.printout = cp2k_properties_linres_polar_print()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&POLAR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.interpolator.status == True:
            self.interpolator.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t&END POLAR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "INTERPOLATOR":
                self.interpolator.set_params({item: params[item]})
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_print_program_run_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_print_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_print_program_run_info_each()

        # basic setting

    def to_input(self, fout):
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

class cp2k_properties_linres_print_restart_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_print_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_print_restart_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END RESTART\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.program_run_info = cp2k_properties_linres_print_program_run_info()
        self.restart = cp2k_properties_linres_print_restart()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.restart.status == True:
            self.restart.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[3] == "RESTART":
                self.restart.set_params({item: params[item]})
            else:
                pass
class cp2k_properties_linres_spinspin_interpolator_conv_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_spinspin_interpolator_conv_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_spinspin_interpolator_conv_info_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&CONV_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END CONV_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_spinspin_interpolator:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.conv_info = cp2k_properties_linres_spinspin_interpolator_conv_info()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&INTERPOLATOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.conv_info.status == True:
            self.conv_info.to_input(fout)
        fout.write("\t\t\t\t&END INTERPOLATOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CONV_INFO":
                self.conv_info.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_spinspin_print_k_matrix_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_spinspin_print_k_matrix:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_spinspin_print_k_matrix_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&POLAR_MATRIX\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END POLAR_MATRIX\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_spinspin_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.polar_matrix = cp2k_properties_linres_spinspin_print_k_matrix()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.polar_matrix.status == True:
            self.polar_matrix.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "K_MATRIX":
                self.polar_matrix.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_spinspin:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.interpolator = cp2k_properties_linres_spinspin_interpolator()
        self.printout = cp2k_properties_linres_spinspin_print()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&SPINSPIN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.interpolator.status == True:
            self.interpolator.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t&END SPINSPIN\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "INTERPOLATOR":
                self.interpolator.set_params({item: params[item]})
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres:
    def __init__(self):
        self.params = {
                "ENERGY_GAP": None,
                "EPS": None,
                "MAX_ITER": None,
                "PRECONDITIONER": None,
                }
        self.status = False

        self.current = cp2k_properties_linres_current()
        self.epr = cp2k_properties_linres_epr()
        self.loclalize = cp2k_properties_linres_localize()
        self.nmr = cp2k_properties_linres_nmr()
        self.polar = cp2k_properties_linres_polar()
        self.printout = cp2k_properties_linres_print()
        self.spinspin = cp2k_properties_linres_spinspin()

        # basic setting
        self.polar.status = True
        self.current.status = True

    def to_input(self, fout):
        fout.write("\t\t&LINRES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.current.status == True:
            self.current.to_input(fout)
        if self.epr.status == True:
            self.epr.to_input(fout)
        if self.localize.status == True:
            self.localize.to_input(fout)
        if self.nmr.status == True:
            self.nmr.to_input(fout)
        if self.polar.status == True:
            self.polar.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.spinspin.status == True:
            self.spinspin.to_input(fout)
        fout.write("\t\t&END LINRES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "CURRENT":
                self.current.set_params({item: params[item]})
            elif item.split("-")[2] == "EPR":
                self.epr.set_params({item: params[item]})
            elif item.split("-")[2] == "LOCALIZE":
                self.localize.set_params({item: params[item]})
            elif item.split("-")[2] == "NMR":
                self.nmr.set_params({item: params[item]})
            elif item.split("-")[2] == "POLAR":
                self.polar.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[2] == "SPINSPIN":
                self.spinspin.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_tddfpt_dipole_moments:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&DIPOLE_MOMENTS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END DIPOLE_MOMENTS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_tddfpt_mgrid:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&MGRID\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END MGRID\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_tddfpt_print_detailed_energy_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
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


class cp2k_properties_tddfpt_print_detailed_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_tddfpt_print_detailed_energy_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&DETAILED_ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END DETAILED_ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_tddfpt_print_guess_vectors_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
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

class cp2k_properties_tddfpt_print_guess_vectors:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_tddfpt_print_guess_vectors_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&GUESS_VECTORS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END GUESS_VECTORS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_tddfpt_print_iteration_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
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


class cp2k_properties_tddfpt_print_iteration_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_tddfpt_print_iteration_info_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&ITERATION_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END ITERATION_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_tddfpt_print_restart_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
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


class cp2k_properties_tddfpt_print_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_tddfpt_print_restart_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END RESTART\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "RESTART":
                self.restart.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_tddfpt_print:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.detailed_energy = cp2k_properties_tddfpt_print_detailed_energy()
        self.guess_vectors = cp2k_properties_tddfpt_print_guess_vectors()
        self.iteration_info = cp2k_properties_tddfpt_print_iteration_info()
        self.restart = cp2k_properties_tddfpt_print_restart()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.detailed_energy.status == True:
            self.detailed_energy.to_input(fout)
        if self.guess_vectors.status == True:
            self.guess_vectors.to_input(fout)
        if self.iteration_info.status == True:
            self.iteration_info.to_input(fout)
        if self.restart.status == True:
            self.restart.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "DETAILED_ENERGY":
                self.detailed_energy.set_params({item: params[item]})
            elif item.split("-")[3] == "GUESS_VECTORS":
                self.guess_vectors.set_params({item: params[item]})
            elif item.split("-")[3] == "ITERATION_INFO":
                self.iteration_info.set_params({item: params[item]})
            elif item.split("-")[3] == "RESTART":
                self.restart.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_tddfpt:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.dipole_moments = cp2k_properties_tddfpt_dipole_moments()
        self.mgrid = cp2k_properties_tddfpt_mgrid()
        self.printout = cp2k_properties_tddfpt_print()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t&TDDFPT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.dipole_moments.status == True:
            self.dipole_moments.to_input(fout)
        if self.mgrid.status == True:
            self.mgrid.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t&END TDDFPT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "DIPOLE_MOMENTS":
                self.dipole_moments.set_params({item: params[item]})
            elif item.split("-")[2] == "MGRID":
                self.mgrid.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass


class cp2k_properties:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.atomic = cp2k_properties_atomic()
        self.et_coupling = cp2k_properties_et_coupling()
        self.fit_charge = cp2k_properties_fit_charge()
        self.resp = cp2k_properties_resp()
        self.linres = cp2k_properties_linres()

    def to_input(self, fout):
        fout.write("\t&PROPERTIES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.atomic.status == True:
            self.atomic.to_input(fout)
        if self.et_coupling.status == True:
            self.et_coupling.to_input(fout)
        if self.fit_charge.status == True:
            self.fit_charge.to_input(fout)
        if self.resp.status == True:
            self.resp.to_input(fout)
        if self.linres.status == True:
            self.linres.to_input(fout)
        fout.write("\t&END PROPERTIES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "ATOMIC":
                self.atomic.set_params({item: params[item]})
            elif item.split("-")[1] == "ET_COUPLING":
                self.et_coupling.set_params({item: params[item]})
            elif item.split("-")[1] == "FIT_CHARGE":
                self.fit_charge.set_params({item: params[item]})
            elif item.split("-")[1] == "RESP":
                self.resp.set_params({item: params[item]})
            elif item.split("-")[1] == "LINRES":
                self.linres.set_params({item: params[item]})
            elif item.split("-")[1] == "TDDFPT":
                self.tddfpt.set_params({item: params[item]})
            else:
                pass
