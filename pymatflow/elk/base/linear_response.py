class casida:
    """
    """
    def __init__(self):
        self.params = {}
        

    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue    
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out
        
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 3:
                self.params[item.split("/")[-1]] = params[item]
                continue
                
class kdotp:
    """
    """
    def __init__(self):
        self.params = {}
        

    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue      
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out
        
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 3:
                self.params[item.split("/")[-1]] = params[item]
                continue
                
class polarizabilities:
    """
    """
    def __init__(self):
        self.params = {}
        

    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue    
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out
        
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 3:
                self.params[item.split("/")[-1]] = params[item]
                continue
                
class scf_in_lr_calculations:
    """
    """
    def __init__(self):
        self.params = {}
        

    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue     
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out
        
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 3:
                self.params[item.split("/")[-1]] = params[item]
                continue
                
class solver:
    """
    """
    def __init__(self):
        self.params = {}
        

    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue     
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out
        
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 3:
                self.params[item.split("/")[-1]] = params[item]
                continue        
        
class static_polarization:
    """
    """
    def __init__(self):
        self.params = {}
        

    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue     
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out
        
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 3:
                self.params[item.split("/")[-1]] = params[item]
                continue        
        
class sternheimer:
    """
    """
    def __init__(self):
        self.params = {}
        

    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue     
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out
        
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 3:
                self.params[item.split("/")[-1]] = params[item]
                continue
                
class vibrational_modes:
    """
    """
    def __init__(self):
        self.params = {}
        
    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue       
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out
        
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 3:
                self.params[item.split("/")[-1]] = params[item]
                continue        


class linear_response:
    """
    """
    def __init__(self):
        self.params = {}
        
        self.casida = casida()
        self.kdotp = kdotp()
        self.polarizabilities = polarizabilities()
        self.scf_in_lr_calculations = scf_in_lr_calculations()
        self.solver = solver()
        self.static_polarization = static_polarization()
        self.sternheimer = sternheimer()
        self.vibrational_modes = vibrational_modes()

    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue    
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        out += self.casida.to_string()
        out += self.kdotp.to_string()
        out += self.polarizabilities.to_string()
        out += self.scf_in_lr_calculations.to_string()
        out += self.solver.to_string()
        out += self.static_polarization.to_string()
        out += self.sternheimer.to_string()
        out += self.vibrational_modes.to_string()
        return out
            
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 2:
                self.params[item.split("/")[-1]] = params[item]
                continue
            if item.split("/")[1] == "Casida":
                self.casida.set_params({item: params[item]})
            elif item.split("/")[1] == "KdotP":
                self.kdotp.set_params({item: params[item]})
            elif item.split("/")[1] == "Polarizabilities":
                self.polarizabilities.set_params({item: params[item]})
            elif item.split("/")[1] == "SCF in LR calculations":
                self.scf_in_lr_calculations.set_params({item: params[item]})
            elif item.split("/")[1] == "Solver":
                self.solver.set_params({item: params[item]})
            elif item.split("/")[1] == "Static Polarization":
                self.static_polarization.set_params({item: params[item]})
            elif item.split("/")[1] == "Sternheimer":
                self.sternheimer.set_params({item: params[item]})
            elif item.split("/")[1] == "Vibrational Modes":
                self.vibrational_modes.set_params({item: params[item]})
            else:
                pass
            
            