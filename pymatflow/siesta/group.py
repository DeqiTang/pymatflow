#from pymatflow.variable import Variable
from pymatflow.variable import VariableGroup
from pymatflow.siesta.base import siesta_variable_to_string

class SiestaVariableGroup(VariableGroup):
    def __init__(self):
        super().__init__()

    #
    def to_string(self):
        out = ""
        for item in self.params:
            out += siesta_variable_to_string(self.params[item])
            out += "\n"
        return out

    def set_params(self, params):
        for key in params:
            self.set_param(key, params[key])