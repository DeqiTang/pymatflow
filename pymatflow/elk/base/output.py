

class berkeleygw:
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

class output:
    """
    """
    def __init__(self):
        self.params = {}
        self.berkeleygw = berkeleygw()

    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue    
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        out += self.berkeleygw.to_string()
        return out
        
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 2:
                self.params[item.split("/")[-1]] = params[item]
                continue
            if item.split("/")[1] == "BerkeleyGW":
                self.berkeleygw.set_params({item: params[item]})
            else:
                pass
            
                            