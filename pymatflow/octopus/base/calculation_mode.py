import sys

class gs:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass


class unocc:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass


class td:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass


class go:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass

class opt_control:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass


class em_resp:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass


class casida:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass

class vdw:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass


class vib_modes:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass


class one_shot:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass


class kdotp:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass


class dummy:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass

class dummy:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass

class invert_ks:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out
    
    def basic_setting(self):
        pass        

class recipe:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass



class calculation_mode:
    """
    """
    def __init__(self):
        self.params = {
            "CalculationMode": None,
        }
        self.gs = gs()
        self.unocc = unocc()
        self.td = td()
        self.go = go()
        self.opt_control =opt_control()
        self.em_resp = em_resp()
        self.casida = casida()
        self.vdw = vdw()
        self.vib_modes = vib_modes()
        self.one_shot = one_shot()
        self.kdotp = kdotp()
        self.dummy = dummy()
        self.invert_ks = invert_ks()
        self.recipe = recipe()

    def to_string(self):
        out  = ""
        #for item in self.params:
        #    out += "%s = %s\n" % (item, self.params[item])
        #    out += "\n"
        #return out
        if self.params["CalculationMode"] == "gs":
            out += self.gs.to_string()
        elif self.params["CalclationMode"] == "unocc":
            out += self.unocc.to_string()
        elif self.params["CalclationMode"] == "td":
            out += self.td.to_string()
        elif self.params["CalclationMode"] == "go":
            out += self.go.to_string()
        elif self.params["CalclationMode"] == "opt_control":
            out += self.opt_control.to_string()
        elif self.params["CalclationMode"] == "em_resp":
            out += self.em_resp.to_string()
        elif self.params["CalclationMode"] == "casida":
            out += self.casida.to_string()
        elif self.params["CalclationMode"] == "vdw":
            out += self.vdw.to_string()
        elif self.params["CalclationMode"] == "vib_modes":
            out += self.vib_modes.to_string()
        elif self.params["CalclationMode"] == "one_shot":
            out += self.one_shot.to_string()
        elif self.params["CalclationMode"] == "kdotp":
            out += self.kdotp.to_string()
        elif self.params["CalclationMode"] == "dummy":
            out += self.dummy.to_string()
        elif self.params["CalclationMode"] == "invert_ks":
            out += self.invert_ks.to_string()
        elif self.params["CalclationMode"] == "recipe":
            out += self.recipe.to_string()


    def basic_setting(self, mode):
        """
        :param: mode-> gs, unocc, td, go, opt_control, em_resp,
                casida, vdw, vib_modes, one_shot, kdotp, dummy,
                invert_ks, recipe
        """
        if mode not in ["gs", "unocc", "td", "go", "opt_control", "em_resp",
                "casida", "vdw", "vib_modes", "one_shot", "kdotp", "dummy",
                "invert_ks", "recipe"]:
            print("===========================================\n")
            print("               WARNING!!!!\n")
            print("octopus.base.calculation_mode.basic_setting:\n")
            print("mode is not in gs, unocc, td, go, opt_control,\n") 
            print("em_resp, casida, vdw, vib_modes, one_shot, kdotp,\n")
            print("dummy, invert_ks, recipe\n")
            sys.exit(1)
        self.params["CalculationMode"] = mode

        if mode == "gs":
            self.gs.basic_setting()
        elif mode == "unocc":
            self.unocc.basic_setting()
        elif mode == "td":
            self.td.basic_setting()
        elif mode == "go":
            self.go.basic_setting()
        elif mode == "opt_control":
            self.opt_control.basic_setting()
        elif mode == "em_resp":
            self.em_resp.basic_setting()
        elif mode == "casida":
            self.casida.basic_setting()
        elif mode == "vdw":
            self.vdw.basic_setting()
        elif mode == "vib_modes":
            self.vib_modes.basic_setting()
        elif mode == "one_shot":
            self.one_shot.basic_setting()
        elif mode == "kdotp":
            self.kdotp.basic_setting()
        elif mode == "dummy":
            self.dummy.basic_setting()
        elif mode == "invert_ks":
            self.invert_ks.basic_setting()
        elif mode == "recipe":
            self.recipe.basic_setting()
