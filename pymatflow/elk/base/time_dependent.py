
class absorbing_boundaries:
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

class photoelectronspectrum:
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

class propagation:
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

class response:
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

class td_output:
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


class time_dependent:
    """
    """
    def __init__(self):
        self.params = {}
        
        self.absorbing_boundaries = absorbing_boundaries()
        self.photoelectronspectrum = photoelectronspectrum()
        self.propagation = propagation()
        self.response = response()
        self.td_output = td_output()

    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue              
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        out += self.absorbing_boundaries.to_string()
        out += self.photoelectronspectrum.to_string()
        out += self.propagation.to_string()
        out += self.response.to_string()
        out += self.td_output.to_string()
        return out
        
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 2:
                self.params[item.split("/")[-1]] = params[item]
                continue
            if item.split("/")[1] == "Absorbing Boundaries":
                self.absorbing_boundaries.set_params({item: params[item]})
            elif item.split("/")[1] == "PhotoElectronSpectrum":
                self.photoelectronspectrum.set_params({item: params[item]})
            elif item.split("/")[1] == "Propagation":
                self.propagation.set_params({item: params[item]})
            elif item.split("/")[1] == "Response":
                self.response.set_params({item: params[item]})
            elif item.split("/")[1] == "TD Output":
                self.td_output.set_params({item: params[item]})
            else:
                pass
                                    