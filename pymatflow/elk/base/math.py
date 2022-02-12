
class rootsolver:
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
                
class sparskit:
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
                
class math:
    """
    """
    def __init__(self):
        self.params = {}
    
        self.rootsolver = rootsolver()
        self.sparskit = sparskit()
    
    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue        
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        out += self.rootsolver.to_string()
        out += self.sparskit.to_string()
        return out
        
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 2:
                self.params[item.split("/")[-1]] = params[item]
                continue
            if item.split("/")[1] == "RootSolver":
                self.rootsolver.set_params({item: params[item]})
            elif item.split("/")[1] == "SPARSKIT":
                self.sparskit.set_params({item: params[item]})
            else:
                pass
            
                    