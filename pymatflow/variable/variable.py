
def n_to_string(n):
    if None == n:
        return ""
    elif n > 0:
        return str(n)
    else:
        return ""


class Variable:
    """
    Note:
        self.value should be always a two dimensional array of string
    """
    def __init__(self, key="", value=""):
        self.status = True
        #self.key = None
        #self.value = []
        self.set(key=key, value=value)
        self.n = 0 # None

    def set(self, key="", value=""):
        self.key = key
        if type(value) is list:
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

    def as_val(self, t=str, shape=0):
        """
        Note:
            return a 2D array of t type if shape == 2
            return a 1D array of t type if shape == 1
            return a scalar of t type if shape == 0
        """
        if 2 == shape:
            out = [[t(val) for val in line] for line in self.value]
        elif 1 == shape:
            out = [t(val) for val in self.value[0]]
        elif 0 == shape:
            out = t(self.value[0][0])
            
    def to_string(self, n=0, indent="", layout="same-line"):
        """
        :param layout:
            'samle-line' or 'second-line'
        """
        if False == self.status:
            return ""
        if "same-line" == layout:
            return self.to_string_same_line(n=n, indent=indent)
        elif "second-line" == layout:
            return self.to_string_second_line(n=n, indent=indent)
        else:
            return ""
        
    def to_string_same_line(self, n=0, indent=""):
        if False == self.status:
            return ""
        
        out = ""
        if 0 == len(self.value):
            return out + self.key
        
        if 1 == len(self.value):
            if 1 == len(self.value[0]):
                out += indent + self.key + n_to_string(n) + " " + self.value[0][0]
            else:
                out += indent + self.key + n_to_string(n)
                for item in self.value[0]:
                    out += " " + item
        else:
            out += indent + self.key + n_to_string(n)
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


    def to_string_second_line(self, n=0, indent=""):
        if False == self.status:
            return ""
        
        out = ""
        
        if 0 == len(self.value):
            return out + self.key

        if 1 == len(self.value):
            if 1 == len(self.value[0]):
                out += indent + self.key + n_to_string(n) + " " + self.value[0][0]
            else:
                out += indent + self.key + n_to_string(n) + "\n"

                out += indent
                for item in self.value[0]:
                    out += " " + item

        else:
            out += indent + self.key + n_to_string(n) + "\n"
            for row in self.value:
                out += indent
                for val in row:
                    out += " " + val
                out += "\n"
        return out


