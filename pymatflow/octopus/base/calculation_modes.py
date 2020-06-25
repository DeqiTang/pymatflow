import sys

class geometry_optimization:
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

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 3:
                self.params[item.split("/")[-1]] = params[item]
                continue


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

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 3:
                self.params[item.split("/")[-1]] = params[item]
                continue

class optimal_control:
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

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 3:
                self.params[item.split("/")[-1]] = params[item]
                continue
                
                

class unoccupied_states:
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

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 3:
                self.params[item.split("/")[-1]] = params[item]
                continue
                

class calculation_modes:
    """
    """
    def __init__(self):
        self.params = {
            "CalculationMode": None,
        }
        self.geometry_optimization = geometry_optimization()
        self.invert_ks = invert_ks()
        self.optimal_contorl = optimal_control()
        self.unoccupied_states = unoccupied_states()

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


    def set_default(self, mode):
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
        pass

    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 2:
                self.params[item.split("/")[-1]] = params[item]
                continue
            if item.split("/")[1] == "Geometry Optimization":
                self.geometry_optimization.set_params({item: params[item]})
            elif item.split("/")[1] == "Invert KS":
                self.invert_ks.set_params({item: params[item]})
            elif item.split("/")[1] == "Optimal Control":
                self.optimal_contorl.set_params({item: params[item]})
            elif item.split("/")[1] == "Unoccupied States":
                slef.unoccupied_states.set_params({item: params[item]})
            
