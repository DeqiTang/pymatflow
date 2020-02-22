#!/usr/bin/env python
# _*_ coding: utf-8 _*_


# ====================
# ====================


class cp2k_dft_efield_constant_env:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&CONSTANT_ENV\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END CONSTANT_ENV\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_dft_efield_custom_env:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&CUSTOM_ENV\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END CUSTOM_ENV\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_efield_gaussian_env:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&GAUSSIAN_ENV\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END GAUSSIAN_ENV\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_efield_ramp_env:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&RAMP_ENV\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RAMP_ENV\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_efield:
    def __init__(self):
        self.params = {
                "ENVELOP": None,
                "INTENSITY": None,
                "PHASE": None,
                "POLARISATION": None,
                "WAVELENGTH": None,
                }
        self.status = False

        self.constant_env = cp2k_dft_efield_constant_env()
        self.custom_env = cp2k_dft_efield_custom_env()
        self.gaussian_env = cp2k_dft_efield_gaussian_env()
        self.ramp_env = cp2k_dft_efield_ramp_env()

    def to_input(self, fout):
        fout.write("\t\t&EFIELD\n")
        for item in self.params:
            if self.params[item] is not None:
                if item == "POLARIZATION":
                    fout.write("\t\t\t%s %.f %.f %.f\n" % (item, self.params[item][0], self.params[item][1], self.params[item][2]))
                else:
                    fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.constant_env.status == True:
            self.constant_env.to_input(fout)
        if self.custom_env.status == True:
            self.cuttom_env.to_input(fout)
        if self.gaussian_env.status == True:
            self.gaussian_env.to_input(fout)
        if self.ramp_env.status == True:
            self.ramp_env.to_input(fout)
        fout.write("\t\t&END EFIELD\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "CONSTANT_ENV":
                self.constant_env.set_params({item: params[item]})
            elif item.split("-")[2] == "CUSTOM_ENV":
                self.custom_env.set_params({item: params[item]})
            elif item.split("-")[2] == "GAUSSIAN_ENV":
                self.gaussian_env.set_params({item: params[item]})
            elif item.split("-")[2] == "RAMP_ENV":
                self.ramp_env.set_params({item: params[item]})
            else:
                pass
