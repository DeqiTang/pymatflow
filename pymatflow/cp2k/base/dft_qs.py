#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ============================
# CP2K / FORCE_EVAL / DFT / QS
# ============================

class cp2k_dft_qs_becke_constraint_atom_group:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&ATOM_GROUP\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END ATOM_GROUP\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_qs_becke_constraint_dummy_atoms:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&DUMMY_ATOMS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END DUMMY_ATOMS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_qs_becke_constraint_program_run_info_each:
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
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_becke_constraint_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_qs_becke_constraint_program_run_info_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END PRQOGRAM_RUN_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_qs_becke_constraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.atom_group = cp2k_dft_qs_becke_constraint_atom_group()
        self.dummy_atoms = cp2k_dft_qs_becke_constraint_dummy_atoms()
        self.program_run_info = cp2k_dft_qs_becke_constraint_program_run_info()

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&BECKE_CONSTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.atom_group.status == True:
            self.atom_group.to_input(fout)
        if self.dummy_atoms.status == True:
            self.dummy.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t\t&END BECKE_CONSTRAINT\n")

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

class cp2k_dft_qs_cdft_hirshfeld_constraint_program_run_info_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_qs_cdft_hirshfeld_constraint_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_qs_cdft_hirshfeld_constraint_program_run_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.to_input(fout)
            else:
                pass

class cp2k_dft_qs_cdft_hirshfeld_constraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.program_run_info = cp2k_dft_qs_cdft_hirshfeld_constraint_program_run_info()
        
        # basic_setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&HIRSHFELD_CONSTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t\t\t&END HIRSHFELD_CONSTRAINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_qs_cdft_outer_scf_cdft_opt:
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
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END CDFT_OPT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_cdft_outer_scf:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cdft_opt = cp2k_dft_qs_cdft_outer_scf_cdft_opt()
        
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&OUTER_SCF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.cdft_opt.status == True:
            self.cdft_opt.to_input(fout)
        fout.write("\t\t\t\t&END OUTER_SCF\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CDFT_OPT":
                self.cdft_opt.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_qs_cdft:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.hirshfeld_constraint = cp2k_dft_qs_cdft_hirshfeld_constraint()
        self.outer_scf = cp2k_dft_qs_cdft_outer_scf()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CDFT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.hirshfeld_constraint.status == True:
            self.hirshfeld_constraint.to_input(fout)
        if self.outer_scf.status == True:
            self.outer_scf.to_input(fout)

        fout.write("\t\t\t&END CDFT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "HIRSHFELD_CONSTRAINT":
                self.hirshfeld_constraint.set_params({item: params[item]})
            elif item.split("-")[3] == "OUTER_SCF":
                self.outer_scf.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_qs_ddapc_restraint_program_run_info_each:
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
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_qs_ddapc_restraint_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.section = "LOW"

        self.each = cp2k_dft_qs_ddapc_restraint_program_run_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.statu == True:
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

class cp2k_dft_qs_ddapc_restraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.program_run_info = cp2k_dft_qs_ddapc_restraint_program_run_info()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&DDAPC_RESTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t\t&END DDAPC_RESTRAINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                if item.split("-")[-1] == "PROGRAM_RUN_INFO":
                    self.program_run_info.section = params[item]
                else:
                    self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_qs_dftb_parameter:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PARAMETER\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END PARAMETER\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_dftb:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.parameter = cp2k_dft_qs_dftb_parameter()

        # basic setting
        self.parameter.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&DFTB\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.parameter.status == True:
            self.parameter.to_input(fout)
        fout.write("\t\t\t&END DFTB\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "PARAMETER":
                self.parameter.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_qs_distribution:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&DISTRIBUTION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END DISTRIBUTION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_lrigpw:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&LRIGPW\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END LRIGPW\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_qs_mulliken_restraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&MULLIKEN_RESTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END MULLIKEN_RESTRAINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_optimize_lri_basis_constraint_exponents:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&CONSTRAINT_EXPONENTS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END CONSTRAINT_EXPONENTS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_optimize_lri_basis:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&OPTIMIZE_LRI_BASIS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END OPTIMIZE_LRI_BASIS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_qs_opt_embed_embed_dens_diff_each:
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
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_qs_opt_embed_embed_dens_diff:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_qs_opt_embed_embed_dens_diff_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EMBED_DENS_DIFF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END EMBED_DENS_DIFF\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.to_input(fout)
            else:
                pass

class cp2k_dft_qs_opt_embed_embed_pot_cube_each:
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
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_opt_embed_embed_pot_cube:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_qs_opt_embed_embed_pot_cube_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EMBED_POT_CUBE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END EMBED_POT_CUBE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.to_input(fout)
            else:
                pass


class cp2k_dft_qs_opt_embed_embed_pot_vector_each:
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
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_qs_opt_embed_embed_pot_vector:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_qs_opt_embed_embed_pot_vector_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EMBED_POT_VECTOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END EMBED_POT_VECTOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_qs_opt_embed:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.embed_dens_diff = cp2k_dft_qs_opt_embed_embed_dens_diff()
        self.embed_pot_cube = cp2k_dft_qs_opt_embed_embed_pot_cube()
        self.embed_pot_vector = cp2k_dft_qs_opt_embed_embed_pot_vector()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&OPT_EMBED\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.embed_dens_diff.status == True:
            self.embed_dens_diff.to_input(fout)
        if self.embed_pot_cube.status == True:
            self.embed_pot_cube.to_input(fout)
        if self.embed_pot_vector.status == True:
            self.embed_pot_vector.to_input(fout)
        fout.write("\t\t\t&END OPT_EMBED\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                if item.split("-")[-1] == "EMBED_DENS_DIFF":
                    self.embed_dens_diff.section = params[item]
                elif item.split("-")[-1] == "EMBED_POT_CUBE":
                    self.embed_pot_cube.section = params[item]
                elif item.split("-")[-1] == "EMBED_POT_VECTOR":
                    self.embed_pot_vector.section = params[item]
                else:
                    self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EMBED_DENS_DIFF":
                self.embed_dens_diff.set_params({item: params[item]})
            elif item.split("-")[3] == "EMBED_POT_CUBE":
                self.embed_pot_cube.set_params({item: params[item]})
            elif item.split("-")[3] == "EMBED_POT_VECTOR":
                self.embed_pot_vector.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_qs_s2_restraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&S2_RESTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END S2_RESTRAINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_se_coulomb:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&COULOMB\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END COULOMB\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_se_exchange:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EXCHANGE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END EXCHANGE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_se_ga:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&GA\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END GA\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_se_lr_correction:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&LR_CORRECTION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END LR_CORRECTION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_qs_se_memory:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&MEMORY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END MEMORY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_qs_se_neighbor_lists:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&NEIGHBOR_LISTS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END NEIGHBOR_LISTS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_se_print_ewald_info_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_se_print_ewald_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_qs_se_print_ewald_info_each()

        # basic seting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EWALD_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END EWALD_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_qs_se_print_neighbor_lists_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_se_print_neighbor_lists:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_qs_se_print_neighbor_lists_each()

        # basic seting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&NEIGHBOR_LISTS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END NEIGHBOR_LISTS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_qs_se_print_subcell_each:
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
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_qs_se_print_subcell:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_qs_se_print_subcell_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&SUBCELL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END SUBCELL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_qs_se_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.ewald_info = cp2k_dft_qs_se_print_ewald_info()
        self.neighbor_lists = cp2k_dft_qs_se_print_neighbor_lists()
        self.subcell = cp2k_dft_qs_se_print_subcell()
    
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.ewald_info.status == True:
            self.ewald_info.to_input(fout)
        if self.neighbor_lists.status == True:
            self.neighbor_lists.to_input(fout)
        if self.subcell.status == True:
            self.subcell.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EWALD_INFO":
                self.ewald_info.set_params({item: params[item]})
            elif item.split("-")[4] == "NEIGHBOR_LISTS":
                self.neighbor_lists.set_params({item: params[item]})
            elif item.split("-")[4] == "SUBCELL":
                self.subcell.set_params({item: params[item]})
            else:
                pass



class cp2k_dft_qs_se_screening:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&SCREENING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END SCREENING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_dft_qs_se:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.coulomb = cp2k_dft_qs_se_coulomb()
        self.exchange = cp2k_dft_qs_se_exchange()
        self.ga = cp2k_dft_qs_se_ga()
        self.lr_correction = cp2k_dft_qs_se_lr_correction()
        self.memory = cp2k_dft_qs_se_memory()
        self.neighbor_lists = cp2k_dft_qs_se_neighbor_lists()
        self.printout = cp2k_dft_qs_se_print()
        self.screening = cp2k_dft_qs_se_screening()

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&SE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.coulomb.status == True:
            self.coulomb.to_input(fout)
        if self.exchange.status == True:
            self.exchange.to_input(fout)
        if self.ga.status == True:
            self.ga.to_input(fout)
        if self.lr_correction.status == True:
            self.lr_correction.to_input(fout)
        if self.memory.status == True:
            self.memory.to_input(fout)
        if self.neighbor_lists.status == True:
            self.neighbor_lists.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.screening.status == True:
            self.screening.to_input(fout)
        fout.write("\t\t\t&END SE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "COULOMB":
                self.coulomb.set_params({item: params[item]})
            elif item.split("-")[3] == "EXCHANGE":
                self.exchange.set_params({item: params[item]})
            elif item.split("-")[3] == "GA":
                self.ga.set_params({item: params[item]})
            elif item.split("-")[3] == "LR_CORRECTION":
                self.lr_correction.set_params({item: params[item]})
            elif item.split("-")[3] == "MEMORY":
                self.memory.set_params({item: params[item]})
            elif item.split("-")[3] == "NEIGHBOR_LISTS":
                self.neighbor_lists.set_params({item: params[item]})
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[3] == "SCREENING":
                self.screening.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_qs:
    def __init__(self):
        self.params = {
                "METHOD": "GPW",
                "EPS_DEFAULT": 1.0E-14,
                "FORCE_PAW": None,
                }
        self.status = False

        self.becke_constraint = cp2k_dft_qs_becke_constraint()
        self.cdft = cp2k_dft_qs_cdft()
        self.ddapc_restraint = cp2k_dft_qs_ddapc_restraint()
        self.dftb = cp2k_dft_qs_dftb()
        self.distribution = cp2k_dft_qs_distribution()
        self.lrigpw = cp2k_dft_qs_lrigpw()
        self.mulliken_restraint = cp2k_dft_qs_mulliken_restraint()
        self.optimize_lri_basis = cp2k_dft_qs_optimize_lri_basis()
        self.opt_embed = cp2k_dft_qs_opt_embed()
        self.s2_restraint = cp2k_dft_qs_s2_restraint()
        self.se = cp2k_dft_qs_se()

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&QS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.becke_constraint.status == True:
            self.becke_constraint.to_input(fout)
        if self.cdft.status == True:
            self.cdft.to_input(fout)
        if self.ddapc_restraint.status == True:
            self.ddapc_restraint.to_input(fout)
        if self.dftb.status == True:
            self.dftb.to_input(fout)
        if self.distribution.status == True:
            self.distribution.to_input(fout)
        if self.lrigpw.status == True:
            self.lrigpw.to_input(fout)
        if self.mulliken_restraint.status == True:
            self.mulliken_restraint.to_input(fout)
        if self.optimize_lri_basis.status == True:
            self.optimize_lri_basis.to_input(fout)
        if self.opt_embed.status == True:
            self.opt_embed.to_input(fout)
        if self.s2_restraint.status == True:
            self.s2_restraint.to_input(fout)
        if self.se.status == True:
            self.se.to_input(fout)
        fout.write("\t\t&END QS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "BECKE_CONSTRAINT":
                self.becke_constraint.set_params({item: params[item]})
            elif item.split("-")[2] == "CDFT":
                self.cdft.set_params({item: params[item]})
            elif item.split("-")[2] == "DDAPC_RESTRAINT":
                sel.fddapc_restraint.set_params({item: params[item]})
            elif item.split("-")[2] == "DFTB":
                self.dftb.set_params({item: params[item]})
            elif item.split("-")[2] == "DISTRIBUTION":
                self.distribution.set_params({item: params[item]})
            elif item.split("-")[2] == "LRIGPW":
                self.lrigpw.set_params({item: params[item]})
            elif item.split("-")[2] == "MULLIKEN_RESTRAINT":
                self.mulliken_restraint.set_params({item: params[item]})
            elif item.split("-")[2] == "OPTIMIZE_LRI_BASIS":
                self.optimize_lri_basis.set_params({item: params[item]})
            elif item.split("-")[2] == "OPT_EMBED":
                self.opt_embed.set_params({item: params[item]})
            elif item.split("-")[2] == "S2_RESTRAINT":
                self.s2_restraint.set_params({item: params[item]})
            elif item.split("-")[2] == "SE":
                self.se.set_params({item: params[item]})
            else:
                pass
                


