#!/usr/bin/env python
# _*_ coding: utf-8 _*_





class cp2k_motion_pint_beads:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&BEADS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END BEADS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_pint_gle:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&GLE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END GLE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_pint_helium:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&HELIUM\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END HELIUM\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
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

class cp2k_motion_pint_nose:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&NOSE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END NOSE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_pint_piglet:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PIGLET\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END PIGLET\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_pint_pile:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PILE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END PILE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_pint_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_pint_qtb:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&QTB\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END QTB\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
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

