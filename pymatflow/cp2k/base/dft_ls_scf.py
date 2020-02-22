#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# =======================================
# CP2K / FORCE_EVAL / DFT /  LS_SCF
# =======================================

class cp2k_dft_ls_scf_chebyshev_dos_each:
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
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&end EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_ls_scf_chebyshev_dos:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_ls_scf_chebyshev_dos_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&DOS\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&end DOS\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_ls_scf_chebyshev_print_specific_e_density_cube_each:
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
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&end EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_ls_scf_chebyshev_print_specific_e_density_cube:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_ls_scf_chebyshev_print_specific_e_density_cube_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PRINT_SPECIFIC_E_DENSITY_CUBE\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&end PRINT_SPECIFIC_E_DENSITY_CUBE\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_ls_scf_chebyshev:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.dos = cp2k_dft_ls_scf_chebyshev_dos()

        self.print_specific_e_density_cube = cp2k_dft_ls_scf_chebyshev_print_specific_e_density_cube()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CHEBYSHEV\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.dos.status == True:
            self.dos.to_input(fout)
        if self.print_specific_e_density_cube.status == True:
            self.print_specific_e_density_cube.to_input(fout)
        fout.write("\t\t\t&end CHEBYSHEV\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "DOS":
                self.dos.set_params({item: params[item]})
            elif item.split("-")[3] == "PRINT_SPECIFIC_E_DENSITY_CUBE":
                self.print_specific_e_density_cube.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_ls_scf_curvy_steps:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CURVY_STEPS\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&end CURVY_STEPS\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_ls_scf_pao_line_search_print_run_info_each:
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
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_ls_scf_pao_line_search_print_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_ls_scf_pao_line_search_print_run_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t&END RUN_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[6] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_ls_scf_pao_line_search_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.run_info = cp2k_dft_ls_scf_pao_line_search_print_run_info()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.run_info.status == True:
            self.run_info.to_input(fout)
        fout.write("\t\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "RUN_INFO":
                self.run_info.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_ls_scf_pao_line_search:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.printout = cp2k_dft_ls_scf_pao_line_search_print()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&LINE_SEARCH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t\t&END LINE_SEARCH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_ls_scf_pao_machine_learning_training_set:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&training_set\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&end training_set\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_ls_scf_pao_machine_learning:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.training_set = cp2k_dft_ls_scf_pao_machine_learning_training_set()
        
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&machine_learning\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.training_set.status == True:
            self.training_set.to_input(fout)
        fout.write("\t\t\t\t&end machine_learning\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "TRAINING_SET":
                self.training_set.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_ls_scf_pao_print_atom_info_each:
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

class cp2k_dft_ls_scf_pao_print_atom_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_ls_scf_pao_print_atom_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&ATOM_INFO\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END ATOM_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass



class cp2k_dft_ls_scf_pao_print_fock_eigenvalues_each:
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

class cp2k_dft_ls_scf_pao_print_fock_eigenvalues:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_ls_scf_pao_print_fock_eigenvalues_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&FOCK_EIGENVALUES\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END FOCK_EIGENVALUES\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_ls_scf_pao_print_fock_gap_each:
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

class cp2k_dft_ls_scf_pao_print_fock_gap:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_ls_scf_pao_print_fock_gap_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&FOCK_GAP\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END FOCK_GAP\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_ls_scf_pao_print_ml_training_data_each:
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

class cp2k_dft_ls_scf_pao_print_ml_training_data:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_ls_scf_pao_print_ml_training_data_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&ML_TRAINING_DATA\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END ML_TRAINING_DATA\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_ls_scf_pao_print_ml_variance_each:
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

class cp2k_dft_ls_scf_pao_print_ml_variance:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_ls_scf_pao_print_ml_variance_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&ML_VARIANCE\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END ML_VARIANCE\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_ls_scf_pao_print_opt_info_each:
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

class cp2k_dft_ls_scf_pao_print_opt_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_ls_scf_pao_print_opt_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&OPT_INFO\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END OPT_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_ls_scf_pao_print_restart_each:
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

class cp2k_dft_ls_scf_pao_print_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_ls_scf_pao_print_restart_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&RESTART\n")
        for item in self.params:
            if self.params[item] is not none:
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
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_ls_scf_pao_print_run_info_each:
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

class cp2k_dft_ls_scf_pao_print_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_ls_scf_pao_print_run_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END RUN_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "each":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_ls_scf_pao_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.atom_info = cp2k_dft_ls_scf_pao_print_atom_info()
        self.fock_eigenvalues = cp2k_dft_ls_scf_pao_print_fock_eigenvalues()
        self.fock_gap = cp2k_dft_ls_scf_pao_print_fock_gap()
        self.ml_training_data = cp2k_dft_ls_scf_pao_print_ml_training_data()
        self.ml_variance = cp2k_dft_ls_scf_pao_print_ml_variance()
        self.opt_info = cp2k_dft_ls_scf_pao_print_opt_info()
        self.restart = cp2k_dft_ls_scf_pao_print_restart()
        self.run_info = cp2k_dft_ls_scf_pao_print_run_info()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.atom_info.status == True:
            self.atom_info.to_input(fout)
        if self.fock_eigenvalues.status == True:
            self.fock_eigenvalues.to_input(fout)
        if self.fock_gap.status == True:
            self.fock_gap.to_input(fout)
        if self.ml_training_data.status == True:
            self.ml_training_data.to_input(fout)
        if self.ml_variance.status == True:
            self.ml_variance.to_input(fout)
        if self.opt_info.status == True:
            self.opt_info.to_input(fout)
        if self.restart.status == True:
            self.restart.to_input(fout)
        if self.run_info.status == True:
            self.run_info.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "ATOM_INFO":
                self.atom_info.set_params({item: params[item]})
            elif item.split("-")[5] == "FOCK_EIGENVALUES":
                self.fock_eigenvalues.set_params({item: params[item]})
            elif item.split("-")[5] == "FOCK_GAP":
                self.fock_gap.set_params({item: params[item]})
            elif item.split("-")[5] == "ML_TRAINING_DATA":
                self.ml_training_data.set_params({item: params[item]})
            elif item.split("-")[5] == "ML_VARIANCE":
                self.ml_variance.set_params({item: params[item]})
            elif item.split("-")[5] == "OPT_INFO":
                self.opt_info.set_params({item: params[item]})
            elif item.split("-")[5] == "RESTART":
                self.restart.set_params({item: params[item]})
            elif item.split("-")[5] == "RUN_INFO":
                self.run_info.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_ls_scf_pao:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.line_search = cp2k_dft_ls_scf_pao_line_search()
        self.machine_learning = cp2k_dft_ls_scf_pao_machine_learning()
        self.printout = cp2k_dft_ls_scf_pao_print()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PAO\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.line_search.status == True:
            self.line_search.to_input(fout)
        if self.machine_learning.status == True:
            self.machine_learning.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t&end PAO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "LINE_SEARCH":
                self.line_search.set_params({item: params[item]})
            elif item.split("-")[3] == "MACHINE_LEARNING":
                self.machine_learning.set_params({item: params[item]})
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_ls_scf_pexsi:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PEXSI\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END PEXSI\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_ls_scf_rho_mixing:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RHO_MIXING\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&end RHO_MIXING\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_dft_ls_scf:
    def __init__(self):
        self.params = {
                "check_s_inv": None,
                "dynamic_threshold": None,
                "eps_diis": None,
                "eps_filter": None,
                "eps_lanczos": None,
                "eps_scf": None,
                "extrapolation_order": None,
                "fixed_mu": None,
                "ini_diis": None,
                "ls_diis": None,
                "matrix_cluster_type": None,
                "max_diis": None,
                "max_iter_lanczos": None,
                "max_scf": None,
                "mixing_fraction": None,
                "mu": None,
                "nmixing": None,
                "non_monotonic": None,
                "perform_mu_scan": None,
                "purification_method": None,
                "report_all_sparsities": None,
                "restart_read": None,
                "restart_write": None,
                "sign_sqrt_order": None,
                "single_precision_matrices": None,
                "s_inversion": None,
                "s_preconditioner": None,
                }
        self.status = False

        self.chebyshev = cp2k_dft_ls_scf_chebyshev()
        self.curvy_steps = cp2k_dft_ls_scf_curvy_steps()
        self.pao = cp2k_dft_ls_scf_pao()
        self.pexsi = cp2k_dft_ls_scf_pexsi()
        self.rho_mixing = cp2k_dft_ls_scf_rho_mixing()

        # basic setting

        self.params["eps_filter"] = 1.0e-7
        self.params["eps_scf"] = 1.0e-5
        self.params["s_preconditioner"] = "atomic"

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&ls_scf\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.chebyshev.status == True:
            self.chebyshev.to_input(fout)
        if self.curvy_steps.status == True:
            self.curvy_steps.to_input(fout)
        if self.pao.status == True:
            self.pao.to_input(fout)
        if self.pexsi.status == True:
            self.pexsi.to_input(fout)
        if self.rho_mixing.status == True:
            self.rho_mixing.to_input(fout)
        fout.write("\t\t&end ls_scf\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "CHEBYSHEV":
                self.chebyshev.set_params({item: params[item]})
            elif item.split("-")[2] == "CURVY_STEPS":
                self.curvy_steps.set_params({item: params[item]})
            elif item.split("-")[2] == "PAO":
                self.pao.set_params({item: params[item]})
            elif item.split("-")[2] == "PEXSI":
                self.pexsi.set_params({item: params[item]})
            elif item.split("-")[2] == "RHO_MIXING":
                self.rho_mixing.set_params({item: params[item]})
            else:
                pass
