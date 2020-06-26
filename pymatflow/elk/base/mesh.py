

class curvilinear:
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
            if self.params[item] == None:
                pass        
            if len(item.split("/")) == 3:
                self.params[item.split("/")[-1]] = params[item]
                continue        
        
class derivatives:
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
                
                
class ffts:
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
                
                
                
class kpoints:
    """
    """
    def __init__(self):
        self.params = {}
        self.kpath = None
        

    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue
            if item == "KPointsGrid":
                out += "%KPointsGrid\n"
                out += "  %d | %d | %d\n" % (self.params[item][0], self.params[item][1], self.params[item][2])
                out += "%\n"
            elif item == "KPointsPath":
                kpath = self.params[item]
                out += "%KPointsPath\n"
                out += "  %d" % kpath[0][4]
                for i in range(1, len(kpath)-1):
                    if kpath[i][4] == "|":
                        out += " | %d" % (0)
                    else:
                        out += " | %d" % (kpath[i][4])
                out += "\n"
                for point in kpath:
                    out += "%f | %f | %f #%s\n" % (point[0], point[1], point[2], point[3])
                out += "%\n"
                out += "\n"
            else:
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
                
                
class simulation_box:
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

class mesh:
    """
    """
    def __init__(self):
        self.params = {}
                
        self.curvilinear = curvilinear()
        self.derivatives = derivatives()
        self.ffts = ffts()
        self.kpoints = kpoints()
        self.simulation_box = simulation_box()
                
        self.kpoints_option = "mp"
        self.kpoints_mp = [1, 1, 1, 0, 0, 0]

    def to_string(self):
        out  = ""
        for item in self.params:
            if self.params[item] == None:
                continue
            if item == "Spacing":
                out += "%Spacing\n"
                out += " %s | %s | %s\n" % (self.params[item][0], self.params[item][1], self.params[item][2])
                out += "%\n"
                out += "\n"
            else:
                out += "%s = %s\n" % (item, self.params[item])
                out += "\n"
        out += self.curvilinear.to_string()
        out += self.derivatives.to_string()
        out += self.ffts.to_string()
        out += self.kpoints.to_string()
        out += self.simulation_box.to_string()
        return out

    def write_kpoints(self, fout):
        """
        :param fout: a file stream for writing
        """
        if self.kpoints_option == "mp":
            fout.write("%KPointsGrid\n")
            fout.write("%d | %d | %d\n" % (
                self.kpoints_mp[0],
                self.kpoints_mp[1],
                self.kpoints_mp[2],
                ))
            fout.write("%\n")                
        elif self.kpoints_option == "kpath":
            # there is a trick:
            # when self.kpath[i][4] == "|"
            # we set the number of k point to connect to the next high symmetry kpoint to 0
            # this is very fantastic !!!
            fout.write("%KPointsPath\n")
            for i in range(len(self.kpath)-2):
                if self.kpath[i][4] == "|":
                    fout.write("0 | ")
                else:
                    fout.write("%d | " % self.kpath[i][4])
                fout.write("%d\n" % self.kpath[i][-2])
            #
            for i in range(len(self.kpath)):
                fout.write("%f | %f | %f #%s\n" % sefl.kpath[i][3])
            #
            #fout.write("KPointsUseSymmetries = no\n")
            fout.write("%\n")      

    def set_kpoints(self, kpoints_mp=[1, 1, 1, 0, 0, 0], option="mp", kpath=None):
        """
        :param kpath: the high symmetry k point path used in bands structure calculation
            in format like this:

            [[kx, ky, kz, label, connect_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', 15], ...]

            if connect_indicator in a kpoint is an integer, then it will connect to the following point
            through the number of kpoints defined by connect_indicator.

            if connect_indicator in a kpoint is '|', then it will not connect to the following point,
        TODO:
        Note:
            "mp" means Monkhorst-Pack scheme
        """
        if option == "mp":
            self.kpoints_option = option
            self.kpoints_mp = kpoints_mp
            return
        if option == "kpath":
            self.kpoints_option = option
            self.kpath = kpath
            return
    
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("/")) == 2:
                self.params[item.split("/")[-1]] = params[item]
                continue
            if item.split("/")[1] == "Curvilinear":
                self.curvilinear.set_params({item: params[item]})
            elif item.split("/")[1] == "Derivatives":
                self.derivatives.set_params({item: params[item]})
            elif item.split("/")[1] == "FFTs":
                self.ffts.set_params({item: params[item]})
            elif item.split("/")[1] == "KPoints":
                self.kpoints.set_params({item: params[item]})
            elif item.split("/")[1] == "Simulation Box":
                self.simulation_box.set_params({item: params[item]})
            else:
                pass            