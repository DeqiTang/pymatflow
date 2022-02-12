from pymatflow.variable import Variable
from pymatflow.variable import VariableGroup
from pymatflow.variable.variable import n_to_string, unit_to_string


def cp2k_variable_to_string_same_line(variable, indent=""):
    """
    Note:
        compared to Variable().to_string_same_line()
        this function will deal with @ which might exist in name of Cp2kSection.name
        or Variable.key. Because in cp2k subsection name and variable can be duplicate.
        And we use %s@%d to distinguish different variable with the same %s key.
        However, when output to string, @%d need to be removed. so this function is defined
        to do that.
    """
    if False == variable.status:
        return ""
    
    if None == variable.value:
        return ""
    
    out = ""
    
    if 0 == len(variable.value):
        return out + variable.key.split("@")[0]
    
    if 1 == len(variable.value):
        if 1 == len(variable.value[0]):
            out += indent + variable.key.split("@")[0] + n_to_string(variable.n) + " " + variable.value[0][0] + " " + unit_to_string(variable.unit)
        else:
            out += indent + variable.key.split("@")[0] + n_to_string(variable.n)
            for item in variable.value[0]:
                out += " " + item
            out += " " + unit_to_string(variable.unit)
    else:
        out += indent + variable.key.split("@")[0] + n_to_string(variable.n) # + " " + unit_to_string(variable.n)
        for val in variable.value[0]:
            out += " " + val
        
        out += "\n"
        for row in range(1, len(variable.value)-1):
            out += indent
            for val in variable.value[row]:
                out += " " + val
            out += "\n"
        out += indent
        for val in variable.value[len(variable.value) - 1]:
            out += " " + val
    return out


class Cp2kSection(VariableGroup):
    def __init__(self, name=""):
        super().__init__()
        self.name = name
        self.section_parameter = ""
        self.section_var = Variable()
        self.status = True
        #self.params = {} # string -> Variable
        self.subsections = {} # string -> Cp2kSection

    def to_string(self, indent="\t"):
        if False == self.status:
            return ""
        out = ""
        
        out += indent + "&" + self.name.split("@")[0] + " " + self.section_parameter + "\n"

        out += cp2k_variable_to_string_same_line(variable=self.section_var, indent=indent+indent) + "\n"

        for item in self.params:
            out += indent + indent + cp2k_variable_to_string_same_line(variable=self.params[item], indent="") + "\n"

        for item in self.subsections:
            out += "\n"
            out += self.subsections[item].to_string(indent+indent)
            out += "\n"

        out += indent + "&end " + self.name.split("@")[0]

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




        





