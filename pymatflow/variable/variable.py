
def n_to_string(n):
    if None == n:
        return ""
    elif n > 0:
        return str(n)
    else:
        return ""

def unit_to_string(unit):
    if None == unit:
        return ""
    else:
        return str(unit)

class Variable:
    """
    Note:
        self.value should be always a two dimensional array of string
    """
    def __init__(self, key="", value="", unit=None):
        self.status = True
        #self.value = []
        self.set(key=key, value=value, unit=unit)
        self.n = 0 # None

    def set(self, key="", value="", unit=None):
        self.key = key
        self.unit = unit
        if None == value:
            self.value = None
            return
        elif type(value) is list:
            if 0 == len(value):
                self.value = [[]]
            else:
                if type(value[0]) is list:
                    self.value = [[str(val) for val in line] for line in value]
                else:
                    self.value = []
                    self.value.append([str(val) for val in value])
        else:
            tmp = [str(value)]
            self.value = []
            self.value.append(tmp)

    def set_n(self, n):
        self.n = n

    def as_val(self, t=str, dim=0):
        """
        Note:
            return a 2D array of t type if dim == 2
            return a 1D array of t type if dim == 1
            return a scalar of t type if dim == 0
        """
        if None == self.value:
            return None

        if 2 == dim:
            out = [[t(val) for val in line] for line in self.value]
        elif 1 == dim:
            out = [t(val) for val in self.value[0]]
        elif 0 == dim:
            out = t(self.value[0][0])
        return out
            
    def to_string(self, indent="", layout="same-line"):
        """
        :param layout:
            'samle-line' or 'second-line'
        """
        if False == self.status:
            return ""
        if "same-line" == layout:
            return self.to_string_same_line(indent=indent)
        elif "second-line" == layout:
            return self.to_string_second_line(indent=indent)
        else:
            return ""
        
    def to_string_same_line(self, indent=""):
        if False == self.status:
            return ""
        
        if None == self.value:
            return ""
        
        out = ""
        
        if 0 == len(self.value):
            return out + self.key
        
        if 1 == len(self.value):
            if 1 == len(self.value[0]):
                out += indent + self.key + n_to_string(self.n) + " " + self.value[0][0] + " " + unit_to_string(self.unit)
            else:
                out += indent + self.key + n_to_string(self.n)
                for item in self.value[0]:
                    out += " " + item
                out += " " + unit_to_string(self.n)
        else:
            out += indent + self.key + n_to_string(self.n) # + " " + unit_to_string(self.n)
            for val in self.value[0]:
                out += " " + val
            
            out += "\n"
            for row in range(1, len(self.value)-1):
                out += indent
                for val in self.value[row]:
                    out += " " + val
                out += "\n"
            out += indent
            for val in self.value[len(self.value) - 1]:
                out += " " + val
        return out


    def to_string_second_line(self, indent=""):
        if False == self.status:
            return ""
        
        if None == self.value:
            return ""
                    
        out = ""
        
        if 0 == len(self.value):
            return out + self.key

        if 1 == len(self.value):
            if 1 == len(self.value[0]):
                out += indent + self.key + n_to_string(self.n) + " " + self.value[0][0] + " " + unit_to_string(self.n)
            else:
                out += indent + self.key + n_to_string(self.n) + "\n"

                out += indent
                for item in self.value[0]:
                    out += " " + item
                OUT += " " + unit_to_string(self.n)
        else:
            out += indent + self.key + n_to_string(self.n) + "\n" + " " + unit_to_string(self.n)
            for row in self.value:
                out += indent
                for val in row:
                    out += " " + val
                out += "\n"
        return out


