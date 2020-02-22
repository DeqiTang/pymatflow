#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# =============================
# CP2K / FORCE_EVAL / DFT / XC
# =============================

class cp2k_dft_xc_adiabatic_rescaling:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&ADIABATIC_RESCALING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END ADIABATIC_RESCALING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_dft_xc_hf_hf_info_each:
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

class cp2k_dft_xc_hf_hf_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xc_hf_hf_info_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&HF_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END HF_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})

class cp2k_dft_xc_hf_interaction_potential:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&INTERACTION_POTENTIAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END INTERACTION_POTENTIAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_hf_load_balance_print_each:
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
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[6] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xc_hf_load_balance_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xc_hf_load_balance_print_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
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

class cp2k_dft_xc_hf_load_balance:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.printout = cp2k_dft_xc_hf_load_balance_print()

        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&LOAD_BALANCE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t\t&END LOAD_BALANCE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_xc_hf_memory:
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
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END MEMORY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_hf_periodic:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PERIODIC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END PERIODIC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_hf_screening:
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
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END SCREENING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_dft_xc_hf:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.hf_info = cp2k_dft_xc_hf_hf_info()
        self.interaction_potential = cp2k_dft_xc_hf_interaction_potential()
        self.load_balance = cp2k_dft_xc_hf_load_balance()
        self.memory = cp2k_dft_xc_hf_memory()
        self.periodic = cp2k_dft_xc_hf_periodic()
        self.screening = cp2k_dft_xc_hf_screening()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&HF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
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
        fout.write("\t\t\t&END HF\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "HF_INFO":
                self.hf_info.set_params({item: parmas[item]})
            elif item.split("-")[3] == "INTERACTION_POTENTIAL":
                self.interaction_potential.set_params({item: params[item]})
            elif item.split("-")[3] == "LOAD_BALANCE":
                self.load_balance.set_params({item: params[item]})
            elif item.split("-")[3] == "MEMORY":
                self.memory.set_params({item: params[item]})
            elif item.split("-")[3] == "PERIODIC":
                self.periodic.set_params({item: params[item]})
            elif item.split("-")[3] == "SCREENING":
                self.screening.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xc_vdw_potential_pair_potential_print_dftd_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_vdw_potential_pair_potential_print_dftd:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xc_vdw_potential_pair_potential_print_dftd_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&PRINT_DFTD\n")
        for item in self.params:
            fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END PRINT_DFTD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xc_vdw_potential_pair_potential:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.print_dftd = cp2k_dft_xc_vdw_potential_pair_potential_print_dftd()

        # basic setting
        self.params["PARAMETER_FILE_NAME"] = "dftd3.dat"

    def to_input(self, fout):
        fout.write("\t\t\t\t&PAIR_POTENTIAL\n")
        for item in self.params:
            fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.print_dftd.status == True:
            self.print_dftd.to_input(fout)
        fout.write("\t\t\t\t&END PAIR_POTENTIAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "PRINT_DFTD":
                self.print_dftd.set_params({item: params[item]})

class cp2k_dft_xc_vdw_potential_non_local:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t&NON_LOCAL\n")
        for item in self.params:
            fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END NON_LOCAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_dft_xc_vdw_potential:
    def __init__(self):
        self.params = {
                "POTENTIAL_TYPE": None,
                }
        self.status = False

        self.pair_potential = cp2k_dft_xc_vdw_potential_pair_potential()
        self.non_local = cp2k_dft_xc_vdw_potential_non_local()
        # basic setting
        self.pair_potential.status = True
        self.params["POTENTIAL_TYPE"] = "PAIR_POTENTIAL"

    def to_input(self, fout):
        fout.write("\t\t\t&VDW_POTENTIAL\n")
        for item in self.params:
            fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.params["POTENTIAL_TYPE"].upper() == "PAIR_POTENTIAL":
            self.pair_potential.to_input(fout)
        elif self.params["POTENTIAL_TYPE"].upper() == "NON_LOCAL":
            self.non_local.to_input(fout)
        fout.write("\t\t\t&END VDW_POTENTIAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "PAIR_POTENTIAL":
                self.pair_potential.set_params({item: params[item]})
            elif item.split("-")[3] == "NON_LOCAL":
                self.non_local.set_params({item: params[item]})



class cp2k_dft_xc_wf_correlation_cphf:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&CPHF\n")
        fout.write("\t\t\t\t&END CPHF\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xc_wf_correlation_direct_canonical:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&WF_DIRECT_CANONICAL\n")
        fout.write("\t\t\t\t&END DIRECT_CANONICAL\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xc_wf_correlation_eri_mme_cutoff_calib:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&CUTOFF_CALIB\n")
        fout.write("\t\t\t\t\t&END CUTOFF_CALIB\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xc_wf_correlation_eri_mme_eri_mme_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xc_wf_correlation_eri_mme_eri_mme_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xc_wf_correlation_eri_mme_eri_mme_info_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&ERI_MME_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END ERI_MME_INFO\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_xc_wf_correlation_eri_mme:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cutoff_calib = cp2k_dft_xc_wf_correlation_eri_mme_cutoff_calib()
        self.eri_mme_info = cp2k_dft_xc_wf_correlation_eri_mme_eri_mme_info()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&ERI_MME\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.cutoff_calib.status == True:
            self.cutoff_calib.to_input(fout)
        if self.eri_mme_info.status == True:
            self.eri_mme_info.to_input(fout)
        fout.write("\t\t\t\t&END ERI_MME\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CUTOFF_CALIB":
                self.cutoff_calib.set_params({item: params[item]})
            elif item.split("-")[4] == "ERI_MME_INFO":
                self.eri_mme_info.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xc_wf_correlation_interaction_potential:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&INTERACTION_POTENTIAL\n")
        fout.write("\t\t\t\t&END INTERACTION_POTENTIAL\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xc_wf_correlation_mp2_info_each:
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
        """
        """
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xc_wf_correlation_mp2_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xc_wf_correlation_mp2_info_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&MP2_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END MP2_INFO\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xc_wf_correlation_opt_ri_basis:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&OPT_RI_BASIS\n")
        fout.write("\t\t\t\t&END OPT_RI_BASIS\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xc_wf_correlation_ri_laplace:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&RI_LAPLACE\n")
        fout.write("\t\t\t\t&END RI_LAPALACE\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xc_wf_correlation_ri_mp2:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&RI_MP2\n")
        fout.write("\t\t\t\t&END RI_MP2\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xc_wf_correlation_ri_rpa_hf_hf_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&EACH\n")
        fout.write("\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xc_wf_correlation_ri_rpa_hf_hf_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xc_wf_correlation_ri_rpa_hf_hf_info_each()

        # basic setting
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&HF_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t&END HF_INFO\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[6] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xc_wf_correlation_ri_rpa_hf_interaction_potential:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&INTERACTION_POTENTIAL\n")
        fout.write("\t\t\t\t\t\t&END INTERACTION_POTENTIAL\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xc_wf_correlation_ri_rpa_hf_load_balance_print_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t\t&EACH\n")
        fout.write("\t\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 9:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xc_wf_correlation_ri_rpa_hf_load_balance_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_xc_wf_correlation_ri_rpa_hf_load_balance_print_each()

        # basic setting 
        self.each.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&LOAD_BALANCE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t&END LOAD_BALANCE\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[6] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xc_wf_correlation_ri_rpa_hf_load_balance:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.printout = cp2k_dft_xc_wf_correlation_ri_rpa_hf_load_balance_print()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&LOAD_BALANCE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t\t\t\t&END LOAD_BALANCE\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[6] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xc_wf_correlation_ri_rpa_hf_memory:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&MEMORY\n")
        fout.write("\t\t\t\t\t\t&END MEMORY\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xc_wf_correlation_ri_rpa_hf_periodic:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&PERIODIC\n")
        fout.write("\t\t\t\t\t\t&END PERIODIC\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xc_wf_correlation_ri_rpa_hf_screening:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&SCREENING\n")
        fout.write("\t\t\t\t\t\t&END SCREENING\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xc_wf_correlation_ri_rpa_hf:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.hf_info = cp2k_dft_xc_wf_correlation_ri_rpa_hf_hf_info()
        self.interaction_potential = cp2k_dft_xc_wf_correlation_ri_rpa_hf_interaction_potential()
        self.load_balance = cp2k_dft_xc_wf_correlation_ri_rpa_hf_load_balance()
        self.memory = cp2k_dft_xc_wf_correlation_ri_rpa_hf_memory()
        self.periodic = cp2k_dft_xc_wf_correlation_ri_rpa_hf_periodic()
        self.screening = cp2k_dft_xc_wf_correlation_ri_rpa_hf_screening()

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&HF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
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
        fout.write("\t\t\t\t\t&end HF\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "HF_INFO":
                self.hf_info.set_params({item: params[item]})
            elif item.split("-")[5] == "INTERACTION_POTENTIAL":
                self.interaction_potential.set_params({item: params[item]})
            elif item.split("-")[5] == "LOAD_BALANCE":
                self.load_balance.set_params({item: params[item]})
            elif item.split("-")[5] == "MEMORY":
                self.memory.set_params({item: params[item]})
            elif item.split("-")[5] == "PERIODIC":
                self.periodic.set_params({item: params[item]})
            elif item.split("-")[5] == "SCREENING":
                self.screening.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_xc_wf_correlation_ri_rpa_im_time:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&IM_TIME\n")
        fout.write("\t\t\t\t\t&end IM_TIME\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_xc_wf_correlation_ri_rpa_ri_axk:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&RI_AXK\n")
        fout.write("\t\t\t\t\t&end RI_AXK\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xc_wf_correlation_ri_rpa_ri_g0w0:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&RI_G0W0\n")
        fout.write("\t\t\t\t\t&end RI_G0W0\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xc_wf_correlation_ri_rpa:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.hf = cp2k_dft_xc_wf_correlation_ri_rpa_hf()
        self.im_time = cp2k_dft_xc_wf_correlation_ri_rpa_im_time()
        self.ri_axk = cp2k_dft_xc_wf_correlation_ri_rpa_ri_axk()
        self.ri_g0w0 = cp2k_dft_xc_wf_correlation_ri_rpa_ri_g0w0()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&RI_RPA\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.hf.status == True:
            self.hf.to_input(fout)
        if self.im_time.status == True:
            self.im_time.to_input(fout)
        if self.ri_axk.status == True:
            self.ri_axk.to_input(fout)
        if self.ri_g0w0.status == True:
            self.ri_g0w0.to_input(fout)
        fout.write("\t\t\t\t&end RI_RPA\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "HF":
                self.hf.set_params({item: params[item]})
            elif item.split("-")[4] == "IM_TIME":
                self.im_time.set_params({item: params[item]})
            elif item.split("-")[4] == "RI_AXK":
                self.ri_axk.set_params({item: params[item]})
            elif item.split("-")[4] == "RI_G0W0":
                self.ri_g0w0.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_xc_wf_correlation_wfc_gpw:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&WFC_GPW\n")
        fout.write("\t\t\t\t&END WFC_GPW\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_xc_wf_correlation:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cphf = cp2k_dft_xc_wf_correlation_cphf()
        self.direct_canonical = cp2k_dft_xc_wf_correlation_direct_canonical()
        self.eri_mme = cp2k_dft_xc_wf_correlation_eri_mme()
        self.interaction_potential = cp2k_dft_xc_wf_correlation_interaction_potential()
        self.mp2_info = cp2k_dft_xc_wf_correlation_mp2_info()
        self.opt_ri_basis = cp2k_dft_xc_wf_correlation_opt_ri_basis()
        self.ri_laplace = cp2k_dft_xc_wf_correlation_ri_laplace()
        self.ri_mp2 = cp2k_dft_xc_wf_correlation_ri_mp2()
        self.ri_rpa = cp2k_dft_xc_wf_correlation_ri_rpa()
        self.wfc_gpw = cp2k_dft_xc_wf_correlation_wfc_gpw()

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&WF_CORRELATION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
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
        fout.write("\t\t\t&END WF_CORRELATION\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "CPHF":
                self.cphf.set_params({item: params[item]})
            elif item.split("-")[3] == "DIRECT_CANONICAL":
                self.direct_canonical.set_params({item: params[item]})
            elif item.split("-")[3] == "ERI_MME":
                self.eri_mme.set_params({item: params[item]})
            elif item.split("-")[3] == "INTERACTION_POTENTIAL":
                self.interaction_potential.set_params({item: params[item]})
            elif item.split("-")[3] == "MP2_INFO":
                self.mp2_info.set_params({item: params[item]})
            elif item.split("-")[3] == "OPT_RI_BASIS":
                self.opt_ri_basis.set_params({item: params[item]})
            elif item.split("-")[3] == "RI_LAPLACE":
                self.ri_laplace.set_params({item: params[item]})
            elif item.split("-")[3] == "RI_MP2":
                self.ri_mp2.set_params({item: params[item]})
            elif item.split("-")[3] == "RI_RPA":
                self.ri_rpa.set_params({item: params[item]})
            elif item.split("-")[3] == "WFC_GPW":
                self.wfc_gpw.set_params({item: params[item]})
            else:
                pass



class cp2k_dft_xc_xc_functional_becke88:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&BECKE88\n")
        fout.write("\t\t\t\t&END BECKE88\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_becke88_lr:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&BECKE88_LR\n")
        fout.write("\t\t\t\t&END BECKE88_LR\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_becke88_lr_adiabatic:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&BECKE88_LR_ADIABATIC\n")
        fout.write("\t\t\t\t&END BECKE88_LR_ADIABATIC\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_dft_xc_xc_functional_becke97:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&BECKE97\n")
        fout.write("\t\t\t\t&END BECKE97\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_becke_roussel:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&BECKE_ROUSSEL\n")
        fout.write("\t\t\t\t&END BECKE_ROUSSEL\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_beef:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&BEEF\n")
        fout.write("\t\t\t\t&END BEEF\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_cs1:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&CS1\n")
        fout.write("\t\t\t\t&END CS1\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_gv09:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&GV09\n")
        fout.write("\t\t\t\t&END GV09\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_hcth:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&HCTH\n")
        fout.write("\t\t\t\t&END HCTH\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_ke_gga:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&KE_GGA\n")
        fout.write("\t\t\t\t&END KE_GGA\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_ke_libxc:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&KE_LIBXC\n")
        fout.write("\t\t\t\t&END KE_LIBXC\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_lda_hole_t_c_lr:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&LDA_HOLE_T_C_LR\n")
        fout.write("\t\t\t\t&END LDA_HOLE_T_C_LR\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_libxc:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&LIBXC\n")
        fout.write("\t\t\t\t&END LIBXC\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_lyp:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&LYP\n")
        fout.write("\t\t\t\t&END LYP\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_lyp_adiabatic:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&LYP_ADIABATIC\n")
        fout.write("\t\t\t\t&END LYP_ADIABATIC\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_optx:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&OPTX\n")
        fout.write("\t\t\t\t&END OPTX\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_p86c:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&P86C\n")
        fout.write("\t\t\t\t&END P86C\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_pade:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PADE\n")
        fout.write("\t\t\t\t&END PADE\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_pbe:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PBE\n")
        fout.write("\t\t\t\t&END PBE\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_pbe_hole_t_c_lr:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PBE_HOLE_T_C_LR\n")
        fout.write("\t\t\t\t&END PBE_HOLE_T_C_LR\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_pw92:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PW92\n")
        fout.write("\t\t\t\t&END PW92\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_pz81:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PZ81\n")
        fout.write("\t\t\t\t&END PZ81\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_tf:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&TF\n")
        fout.write("\t\t\t\t&END TF\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_tfw:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&TFW\n")
        fout.write("\t\t\t\t&END TFW\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_tpss:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&TPSS\n")
        fout.write("\t\t\t\t&END TPSS\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_vwn:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&VWN\n")
        fout.write("\t\t\t\t&END VWN\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_xalpha:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&XALPHA\n")
        fout.write("\t\t\t\t&END XALPHA\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_xgga:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&XGGA\n")
        fout.write("\t\t\t\t&END XGGA\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional_xwpbe:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&XWPBE\n")
        fout.write("\t\t\t\t&END XWPBE\n")

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_functional:
    def __init__(self):
        self.section = "PBE"
        self.params = {
                }
        self.status = False

        self.becke88 =  cp2k_dft_xc_xc_functional_becke88()
        self.becke88_lr = cp2k_dft_xc_xc_functional_becke88_lr()
        self.becke88_lr_adiabatic = cp2k_dft_xc_xc_functional_becke88_lr_adiabatic()
        self.becke97 = cp2k_dft_xc_xc_functional_becke97()
        self.becke_roussel = cp2k_dft_xc_xc_functional_becke_roussel()
        self.beef = cp2k_dft_xc_xc_functional_beef()
        self.cs1 = cp2k_dft_xc_xc_functional_cs1()
        self.gv09 = cp2k_dft_xc_xc_functional_gv09()
        self.hcth = cp2k_dft_xc_xc_functional_hcth()
        self.ke_gga = cp2k_dft_xc_xc_functional_ke_gga()
        self.ke_libxc = cp2k_dft_xc_xc_functional_ke_libxc()
        self.lda_hole_t_c_lr = cp2k_dft_xc_xc_functional_lda_hole_t_c_lr()
        self.libxc = cp2k_dft_xc_xc_functional_libxc()
        self.lyp = cp2k_dft_xc_xc_functional_lyp()
        self.lyp_adiabatic = cp2k_dft_xc_xc_functional_lyp_adiabatic()
        self.optx = cp2k_dft_xc_xc_functional_optx()
        self.p86c = cp2k_dft_xc_xc_functional_p86c()
        self.pade = cp2k_dft_xc_xc_functional_pade()
        self.pbe = cp2k_dft_xc_xc_functional_pbe()
        self.pbe_hole_t_c_lr = cp2k_dft_xc_xc_functional_pbe_hole_t_c_lr()
        self.pw92 = cp2k_dft_xc_xc_functional_pw92()
        self.pz81 = cp2k_dft_xc_xc_functional_pz81()
        self.tf = cp2k_dft_xc_xc_functional_tf()
        self.tfw = cp2k_dft_xc_xc_functional_tfw()
        self.tpss = cp2k_dft_xc_xc_functional_tpss()
        self.vwn = cp2k_dft_xc_xc_functional_vwn()
        self.xalpha = cp2k_dft_xc_xc_functional_xalpha()
        self.xgga = cp2k_dft_xc_xc_functional_xgga()
        self.xwpbe = cp2k_dft_xc_xc_functional_xwpbe()

        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&XC_FUNCTIONAL %s\n" % self.section)
        
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

        fout.write("\t\t\t&END XC_FUNCTIONAL\n")

    def set_params(self, params):
        """
        set_params for xc_functional is different from many other
        set_params, as it deal with the key 'XC_FUNCTIONAL' only
        """
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "BECKE88":
                self.becke88.set_params({item: params[item]})
            elif item.split("-")[3] == "BECKE88_LR":
                self.becke88_lr.set_params({item: params[item]})
            elif item.split("-")[3] == "BECKE88_LR_ADIABATIC":
                self.becke88_lr_adiabatic.set_params({item: params[item]})
            elif item.split("-")[3] == "BECKE97":
                self.becke97.set_params({item: params[item]})
            elif item.split("-")[3] == "BECKE_ROUSSEL":
                self.becke_roussel.set_params({item: params[item]})
            elif item.split("-")[3] == "BEEF":
                self.beef.set_params({item: params[item]})
            elif item.split("-")[3] == "CS1":
                self.cs1.set_params({item: params[item]})
            elif item.split("-")[3] == "GV09":
                self.gv09.set_params({item: params[item]})
            elif item.split("-")[3] == "HCTH":
                self.hcth.set_params({item: params[item]})
            elif item.split("-")[3] == "KE_GGA":
                self.ke_gga.set_params({item: params[item]})
            elif item.split("-")[3] == "KE_LIBXC":
                self.ke_libxc.set_params({item: params[item]})
            elif item.split("-")[3] == "LDA_HOLE_T_C_LR":
                self.lda_hole_t_c_lr.set_params({item: params[item]})
            elif item.split("-")[3] == "LIBXC":
                self.libxc.set_params({item: params[item]})
            elif item.split("-")[3] == "LYP":
                self.lyp.set_params({item: params[item]})
            elif item.split("-")[3] == "LYP_ADIABATIC":
                self.lyp_adiabatic.set_params({item: params[item]})
            elif item.split("-")[3] == "OPTX":
                self.optx.set_params({item: params[item]})
            elif item.split("-")[3] == "P86C":
                self.p86c.set_params({item: params[item]})
            elif item.split("-")[3] == "PADE":
                self.pade.set_params({item: params[item]})
            elif item.split("-")[3] == "PBE":
                self.pbe.set_params({item: params[item]})
            elif item.split("-")[3] == "PBE_HOLE_T_C_LR":
                self.pbe_hole_t_c_lr.set_params({item: params[item]})
            elif item.split("-")[3] == "PW92":
                self.pw92.set_params({item: params[item]})
            elif item.split("-")[3] == "PZ81":
                self.pz81.set_params({item: params[item]})
            elif item.split("-")[3] == "TF":
                self.tf.set_params({item: params[item]})
            elif item.split("-")[3] == "TFW":
                self.tfw.set_params({item: params[item]})
            elif item.split("-")[3] == "TPSS":
                self.tpss.set_params({item: params[item]})
            elif item.split("-")[3] == "VWN":
                self.vwn.set_params({item: params[item]})
            elif item.split("-")[3] == "XALPHA":
                self.xalpha.set_params({item: params[item]})
            elif item.split("-")[3] == "XGGA":
                self.xgga.set_params({item: params[item]})
            elif item.split("-")[3] == "XWPBE":
                self.xwpbe.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_xc_xc_grid:
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
        fout.write("\t\t\t&XC_GRID\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END XC_GRID\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]


class cp2k_dft_xc_xc_potential_saop:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&SAOP\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END SAOP\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]

class cp2k_dft_xc_xc_potential:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.saop = cp2k_dft_xc_xc_potential_saop()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&XC_POTENTIAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.saop.status == True:
            self.saop.to_input(fout)
        fout.write("\t\t\t&END XC_POTENTIAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "SAOP":
                self.saop.set_params({item: params[item]})


class cp2k_dft_xc:
    def __init__(self):
        self.params = {
                "DENSITY_CUTOFF": None,
                "DENSITY_SMOOTH_CUTOFF_RANGE": None,
                "FUNCTIONAL_ROUTINE": None,
                "GRADIENT_CUTOFF": None,
                "TAU_CUTOFF": None,
                }
        self.status = False

        self.adiabatic_rescaling = cp2k_dft_xc_adiabatic_rescaling()
        self.hf = cp2k_dft_xc_hf()
        self.vdw_potential = cp2k_dft_xc_vdw_potential()
        self.wf_correlation = cp2k_dft_xc_wf_correlation()
        self.xc_functional = cp2k_dft_xc_xc_functional()
        self.xc_grid = cp2k_dft_xc_xc_grid()
        self.xc_potential = cp2k_dft_xc_xc_potential()

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
        fout.write("\t\t&XC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s" % (item, str(self.params[item])))
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
        fout.write("\t\t&END XC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                if item.split("-")[-1] == "XC_FUNCTIONAL":
                    self.xc_functional.section = params[item]
                elif item.split("-")[-1] == "VDW_POTENTIAL":
                    self.vdw_potential.section = params[item]
                else:
                    self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "XC_FUNCTIONAL" and len(item.split("-")) > 3:
                self.xc_functional.set_params({item: params[item]})
            elif item.split("-")[2] == "VDW_POTENTIAL":
                self.vdw_potential.set_params({item: params[item]})
            elif item.split("-")[2] == "ADIABATIC_RESCALING":
                self.adiabatic_rescaling.set_params({item: params[item]})
            elif item.split("-")[2] == "HF":
                self.hf.set_params({item: params[item]})
            elif item.split("-")[2] == "WF_CORRELATION":
                self.wf_correlation.set_params({item: params[item]})
            elif item.split("-")[2] == "XC_GRID":
                self.xc_grid.set_params({item: params[item]})
            elif item.split("-")[2] == "XC_POTENTIAL":
                self.xc_potential.set_params({item: params[item]})
            else:
                pass

