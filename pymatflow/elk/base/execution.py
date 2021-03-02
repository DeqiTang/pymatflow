
import sys
import os

class accel:
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
                
class debug:
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
                                
class io:
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
                                                
class optimization:
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
                                                                
class parallelization:
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
                                                                                
class symmetries:
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
                                                                                                
class units:
    def __init__(self):
        self.params = {}
        
        self.params["Units"] = "ev_angstrom"
        self.params["UnitsOutput"] = "ev_angstrom"

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
                                                                                                                

class execution:
    """
    """
    def __init__(self):
        self.params = {}
        
        self.accel = accel()
        self.debug = debug()
        self.io = io()
        self.optimization = optimization()
        self.parallelization = parallelization()
        self.symmetries = symmetries()
        self.units = units()

    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue        
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        out += self.accel.to_string()
        out += self.debug.to_string()
        out += self.io.to_string()
        out += self.optimization.to_string()
        out += self.parallelization.to_string()
        out += self.symmetries.to_string()
        out += self.units.to_string()
        return out
        
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 2:
                self.params[item.split("/")[-1]] = params[item]
                continue
            if item.split("/")[1] == "Accel":
                self.accel.set_params({item: params[item]})
            elif item.split("/")[1] == "Debug":
                self.debug.set_params({item: params[item]})
            elif item.split("/")[1] == "IO":
                self.io.set_params({item: params[item]})
            elif item.split("/")[1] == "Optimization":
                self.optimization.set_params({item: params[item]})
            elif item.split("/")[1] == "Parallelization":
                self.parallelization.set_params({item: params[item]})
            elif item.split("/")[1] == "Symmetries":
                self.symmetries.set_params({item: params[item]})
            elif item.split("/")[1] == "Units":
                slef.units.set_params({item: params[item]})
            
