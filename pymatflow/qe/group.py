#from pymatflow.variable import Variable
from pymatflow.variable import VariableGroup

class QeVariableGroup(VariableGroup):
    def __init__(self):
        super().__init__()

    #
    def to_string(self):
        out = ""
        for item in self.params:
            out += self.params[item].to_string()
            out += "\n"
        return out
        
    def set_params(self, params):
        for key in params:
            self.set_param(key, params[key])