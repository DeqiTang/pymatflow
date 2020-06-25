
class oct_casida_spectrum:
    """
    """
    def __init__(self):
        self.params = {}
        
    def to_string(self):
        out  = ""
        for item in self.params:
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


class oct_center_geom:
    """
    """
    def __init__(self):
        self.params = {}
        
    def to_string(self):
        out  = ""
        for item in self.params:
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


class oct_conductivity_spectrum:
    """
    """
    def __init__(self):
        self.params = {}
        
    def to_string(self):
        out  = ""
        for item in self.params:
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


class oct_convert:
    """
    """
    def __init__(self):
        self.params = {}
        
    def to_string(self):
        out  = ""
        for item in self.params:
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


class oct_local_multipoles:
    """
    """
    def __init__(self):
        self.params = {}
        
    def to_string(self):
        out  = ""
        for item in self.params:
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


class oct_photoelectron_spectrum:
    """
    """
    def __init__(self):
        self.params = {}
        
    def to_string(self):
        out  = ""
        for item in self.params:
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


class oct_propagation_spectrum:
    """
    """
    def __init__(self):
        self.params = {}
        
    def to_string(self):
        out  = ""
        for item in self.params:
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


class oct_test:
    """
    """
    def __init__(self):
        self.params = {}
        
    def to_string(self):
        out  = ""
        for item in self.params:
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


class oct_vibrational_spectrum:
    """
    """
    def __init__(self):
        self.params = {}
        
    def to_string(self):
        out  = ""
        for item in self.params:
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


class oct_xyz_anim:
    """
    """
    def __init__(self):
        self.params = {}
        
    def to_string(self):
        out  = ""
        for item in self.params:
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



class utilities:
    """
    """
    def __init__(self):
        self.params = {}

        self.oct_casida_spectrum = oct_casida_spectrum()
        self.oct_center_geom = oct_center_geom()
        self.oct_conductivity_spectrum = oct_conductivity_spectrum()
        self.oct_convert = oct_convert()
        self.oct_local_mulitipoles = oct_local_multipoles()
        self.oct_photoelectron_spectrum = oct_photoelectron_spectrum()
        self.oct_propagation_spectrum = oct_propagation_spectrum()
        self.oct_test = oct_test()
        self.oct_vibrational_spectrum = oct_vibrational_spectrum()
        self.oct_xyz_anim = oct_xyz_anim()

    def to_string(self):
        out  = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out
        
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 2:
                self.params[item.split("/")[-1]] = params[item]
                continue
            if item.split("/")[1] == "oct-casida_spectrum":
                self.oct_casida_spectrum.set_params({item: params[item]})
            elif item.split("/")[1] == "oct-center-geom":
                self.oct_center_geom.set_params({item: params[item]})
            elif item.split("/")[1] == "oct-conductivity_spectrum":
                self.oct_conductivity_spectrum.set_params({item: params[item]})
            elif item.split("/")[1] == "oct-convert":
                self.oct_convert.set_params({item: params[item]})
            elif item.split("/")[1] == "oct-local_multipoles":
                self.oct_local_multipoles.set_params({item: params[item]})
            elif item.split("/")[1] == "oct-photoelectron_spectrum":
                self.oct_photoelectron_spectrum.set_params({item: params[item]})
            elif item.split("/")[1] == "oct-propagation_spectrum":
                self.oct_propagation_spectrum.set_params({item: params[item]})
            elif item.split("/")[1] == "oct-test":
                self.oct_test.set_params({item: params[item]})
            elif item.split("/")[1] == "oct-vibrational_spectrum":
                self.oct_vibrationl_spectrm.set_params({item: params[item]})
            elif item.split("/")[1] == "oct-xyz-anim":
                self.oct_xyz_anim.set_params({item: params[item]})
            else:
                pass
                                          