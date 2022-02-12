import sys

class geometry_optimization:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            if self.params[item] == None:
                continue      
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
            if self.params[item] == None:
                continue       
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
            if self.params[item] == None:
                continue       
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def basic_setting(self):
        pass

    def set_params(self, params):
        """
        """
        for item in params:
            if self.params[item] == None:
                pass        
            if len(item.split("/")) == 3:
                self.params[item.split("/")[-1]] = params[item]
                continue
                
                

class unoccupied_states:
    def __init__(self):
        self.params = {}

    def to_string(self):
        out = ""
        for item in self.params:
            if self.params[item] == None:
                continue      
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
        self.optimal_control = optimal_control()
        self.unoccupied_states = unoccupied_states()

    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue            
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        
        out += self.geometry_optimization.to_string()
        out += self.invert_ks.to_string()
        out += self.optimal_control.to_string()
        out += self.unoccupied_states.to_string()
        return out

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
                self.optimal_control.set_params({item: params[item]})
            elif item.split("/")[1] == "Unoccupied States":
                slef.unoccupied_states.set_params({item: params[item]})
            
