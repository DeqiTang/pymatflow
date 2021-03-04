from pymatflow.abinit.group import AbinitVariableGroup
"""
in control of ions moving related parameters
"""

class AbinitIons(AbinitVariableGroup):
    """
    """
    def __init__(self):
        super().__init__()
        self.incharge = [
            "ionmov", "optcell", "ecutsm", "ntime", "tolmxde", "tolmxf",
            "dilatmx",
            ]
        self.status = True
    def to_input(self, fout):
        """
        :param fout: a file stream for writing
        """
        fout.write(self.to_string())
        #

    def to_string(self, n=0):
        """
        :return input_str is the string of all the set params
        """
        input_str = ""
        input_str += "# ============================\n"
        input_str += "# ions moving related settting\n"
        input_str += "# ============================\n"
        input_str += "\n"
        self.set_n(n)
        input_str += super().to_string()
        return input_str
        #

    def basic_setting(self, mode="opt"):
        """
        :param mode: opt or md
        """
        if mode == "opt":
            self.set_param("ionmov", 3)
            self.set_param("optcell", 0)
            self.set_param("ntime", 100)
            #self.set_param("nctime", 1) # for netcdf writing
            self.set_param("tolmxf", 0)  #5.0e-4 # Ha/Bohr
            self.set_param("tolmxde", "0.0001 eV")
        elif mode == "md":
            self.set_param("ionmov", 8) #11
            self.set_param("dtion", 100)
            self.set_param("ntime", 1000)
            self.set_param("nctime", 1)
            self.set_param("mdtemp(1)", 300)
            self.set_param("mdtemp(2)", 300)
            self.set_param("tolmxf", 5.0e-4) # Ha/Bohr


class AbinitIonsOld:
    """
    """
    def __init__(self):
        self.params = {}
        self.incharge = [
            "ionmov", "optcell", "ecutsm", "ntime", "tolmxde", "tolmxf",
            "dilatmx",
            ]
        self.status = True
    def to_input(self, fout):
        """
        :param fout: a file stream for writing
        """
        fout.write("# ============================\n")
        fout.write("# ions moving related settting\n")
        fout.write("# ============================\n")
        fout.write("\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item ,str(self.params[item])))
                fout.write("\n")
        fout.write("\n\n")
        #

    def to_string(self, n=0):
        """
        :return input_str is the string of all the set params
        """
        input_str = ""
        input_str += "# ============================\n"
        input_str += "# ions moving related settting\n"
        input_str += "# ============================\n"
        input_str += "\n"
        for item in self.params:
            if self.params[item] is not None:
                input_str += "%s%s %s\n" % (item, n if n > 0 else "", str(self.params[item]))
                input_str += "\n"
        input_str += "\n\n"
        return input_str
        #

    def basic_setting(self, mode="opt"):
        """
        :param mode: opt or md
        """
        if mode == "opt":
            self.params["ionmov"] = 3
            self.params["optcell"] = 0
            self.params["ntime"] = 100
            #self.params["nctime"] = 1 # for netcdf writing
            self.params["tolmxf"] = 0 #5.0e-4 # Ha/Bohr
            self.params["tolmxde"] = "0.0001 eV"
        elif mode == "md":
            self.params["ionmov"] = 8 #11
            self.params["dtion"] = 100
            self.params["ntime"] = 1000
            self.params["nctime"] = 1
            self.params["mdtemp(1)"] = 300
            self.params["mdtemp(2)"] = 300
            self.params["tolmxf"] = 5.0e-4 # Ha/Bohr

    def set_params(self, params):
        for item in params:
            self.params[item] = params[item]
