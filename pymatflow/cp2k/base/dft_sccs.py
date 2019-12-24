#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ==================
# ==================

class cp2k_dft_sccs_andreussi:
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&ANDREUSSI\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END ANDREUSSI\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_sccs_fattebert_gygi:
    """
    Note:
        FATTEBERT-GYGI is FATTEBERT-GYGI insread of FATTEBERT_GYGI
        but when using set_params to set it we still use:
        FATTTEBERT_GYGI to keep consistency.
    """
    def __init__(self):
        self.params = {}
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&FATTEBERT-GYGI\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END FATTEBERT-GYGI\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_dft_sccs:
    def __init__(self):
        self.params = {}
        self.status = False

        self.andreussi = cp2k_dft_sccs_andreussi()
        self.fattebert_gygi = cp2k_dft_sccs_fattebert_gygi()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&SCCS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.andreussi.status == True:
            self.andreussi.to_input(fout)
        if self.fattebert_gygi.status == True:
            self.fattebert_gygi.to_input(fout)
        fout.write("\t\t&END SCCS\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "ANDREUSSI":
                self.andreussi.set_params({item: params[item]})
            elif item.split("-")[2] == "FATTEBERT_GYGI":
                self.fattebert_gygi.set_params({item: params[item]})
            else:
                pass


