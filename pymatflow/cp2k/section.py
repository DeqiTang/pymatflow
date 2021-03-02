from pymatflow.variable import Variable


class Cp2kSection:
    def __init__(self, name=""):
        self.name = name
        self.section_parameter = ""
        self.section_var = Variable()
        self.status = True
        self.params = {} # string -> Variable
        self.subsections = {} # string -> Cp2kSection

    def to_string(self, indent="\t"):
        if False == self.status:
            return ""
        out = ""
        
        out += indent + "&" + self.name + " " + self.section_parameter + "\n"

        out += self.section_var.to_string(layout="same-line", indent=indent+indent) + "\n"

        for item in self.params:
            out += indent + indent + self.params[item].to_string(layout="same-line", indent="") + "\n"

        for item in self.subsections:
            out += "\n"
            out += self.subsections[item].to_string(indent+indent)
            out += "\n"

        out += indent + "&end " + self.name
        return out


    def add_subsection(self, key, section=None):
        if None == section:
            self.subsections[key] = Cp2kSection(key)
        else:
            self.subsections[key] = section
            self.subsections[key].name = key
        return self.subsections[key]

    def remove_subsection(self, key):
        if key in self.subsections:
            del self.subsections[key]

    def set_param(self, key, value):
        self.remove(key)
        self.params[key] = Variable(key, value)

    def contains(self, key):
        if key in self.params:
            return True
        else:
            return False

    def set_status(self, key, status):
        if False == self.contains(key):
            return
        else:
            self.params[key].status = status

    def remove(self, key):
        if key in self.params:
            del self.params[key]

    def clear(self):
        self.params.clear()

    def get(self, key, t=str, shape=0):
        """
        Note:
            return a 2D array of t type if shape == 2
            return a 1D array of t type if shape == 1
            return a scalar of t type if shape == 0
        """
        return self.params[key].as_val(t=t, shape=shape)


        





