#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ==================================
# CP2K / FORCE_EVAL / DFT / ALMO_SCF
# ==================================

class cp2k_dft_almo_scf_almo_optimizer_diis:
    def __init__(self):
        self.params = {

                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&ALMO_OPTIMIZAER_DIIS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END ALMO_OPTIMIZER_DIIS\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_almo_scf_almo_optimizer_pcg:
    def __init__(self):
        self.params = {

                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&ALMO_OPTIMIZAER_PCG\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END ALMO_OPTIMIZER_PCG\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_almo_scf_analysis_print_almo_cta_each:
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


class cp2k_dft_almo_scf_analysis_print_almo_cta:
    def __init__(self):
        self.params = {

                }
        self.status = False

        self.each = cp2k_dft_almo_scf_analysis_print_almo_cta_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&ALMO_CTA\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END ALMO_CTA\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_almo_scf_analysis_print_almo_eda_ct_each:
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

class cp2k_dft_almo_scf_analysis_print_almo_eda_ct:
    def __init__(self):
        self.params = {

                }
        self.status = False

        self.each = cp2k_dft_almo_scf_analysis_print_almo_eda_ct_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&ALMO_EDA_CT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END ALMO_EDA_CT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_almo_scf_analysis_print:
    def __init__(self):
        self.params = {

                }
        self.status = False

        self.almo_cta = cp2k_dft_almo_scf_analysis_print_almo_cta()
        self.almo_eda_ct = cp2k_dft_almo_scf_analysis_print_almo_eda_ct()

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.almo_cta.status == True:
            self.almo_cat.to_input(fout)
        if self.almo_eda_ct.status == True:
            self.almo_eda_ct.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "ALMO_CTA":
                self.almo_cta.set_params({item: params[item]})
            elif item.split("-")[4] == "ALMO_EDA_CT":
                self.almo_eda_ct.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_almo_scf_analysis:
    def __init__(self):
        self.params = {

                }
        self.status = False

        self.printout = cp2k_dft_almo_scf_analysis_print()
        
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&ANALYSIS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t&END ANALYSIS\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_almo_scf_xalmo_optimizer_pcg:
    def __init__(self):
        self.params = {

                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&XALMO_OPTIMIZER_PCG\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END XALMO_OPTIMIZER_PCG\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_almo_scf:
    def __init__(self):
        self.params = {

                }
        self.status = False

        self.almo_optimizer_diis = cp2k_dft_almo_scf_almo_optimizer_diis()
        self.almo_optimizer_pcg = cp2k_dft_almo_scf_almo_optimizer_pcg()
        self.anlysis = cp2k_dft_almo_scf_analysis()
        self.xalmo_optimizer_pcg = cp2k_dft_almo_scf_xalmo_optimizer_pcg()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&ALMO_SCF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.almo_optimizer_diis.status == True:
            self.almo_optimizer_diis.to_input(fout)
        if self.almo_optimizer_pcg.status == True:
            self.almo_optimizer_pcg.to_input(fout)
        if self.analysis.status == True:
            self.analysis.to_input(fout)
        if self.xalmo_optimizer_pcg.status == True:
            self.xalmo_optimizer_pcg.to_input(fout) 
        fout.write("\t\t&END ALMO_SCF\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "ALMO_OPTIMIZER_DIIS":
                self.almo_optimizer_diis.set_params({item: params[item]})
            elif item.split("-")[2] == "ALMO_OPTIMIZER_PCG":
                self.almo_optimizer_pcg.set_params({item: params[item]})
            elif item.split("-")[2] == "ANLYSIS":
                self.analysis.set_params({item: params[item]})
            elif item.split("-")[2] == "XALMO_OPTIMIZER_PCG":
                self.xalmo_optimizer_pcg.set_params({item: params[item]})
            else:
                pass
