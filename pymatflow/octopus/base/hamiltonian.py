
class dft_u:
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

class pcm:
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

class poisson:
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

class xc:
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


class hamiltonian:
    """
    """
    def __init__(self):
        self.params = {}

        self.dft_u = dft_u()
        self.pcm = pcm()
        self.poisson = poisson()
        self.xc = xc()
        
    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue       
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        out += self.dft_u.to_string()
        out += self.pcm.to_string()
        out += self.poisson.to_string()
        out += self.xc.to_string()
        return out
        
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 2:
                self.params[item.split("/")[-1]] = params[item]
                continue
            if item.split("/")[1] == "DFT+U":
                self.dft_u.set_params({item: params[item]})
            elif item.split("/")[1] == "PCM":
                self.pcm.set_params({item: params[item]})
            elif item.split("/")[1] == "Poisson":
                self.poisson.set_params({item: params[item]})
            elif item.split("/")[1] == "XC":
                self.xc.set_params({item: params[item]})
            else:
                pass
            
                    