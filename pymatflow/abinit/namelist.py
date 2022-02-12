#from pymatflow.variable import Variable
from pymatflow.variable import VariableGroup

class OpticNamelist(VariableGroup):
    def __init__(self, name=None):
        super().__init__()
        self.name = name

    def to_string(self):
        out = "&%s\n" % self.name
        for item in self.params:
            if len(self.params[item].value) == 1:
                if len(self.params[item].value[0]) == 1:
                    if item in ["ddkfile_1", "ddkfile_2", "ddkfile_3", "wfkfile"]:
                        out += "%s = \'%s\'" % (item, self.params[item].value[0][0])
                    else:
                        out += "%s = %s" % (item, self.params[item].as_val(t=str))
                else:
                    out += "%s = %s" %  (item, str(self.params[item].value[0][0]))
                    for other in self.params[item].value[0][1:]:
                        out += ", %s" % str(other)
            else:
                out += "%s -> 2D matrix is not supported now" % item
            out += ","
            out += "\n"
        out += "/\n"

        return out
        
    def set_params(self, params):
        for key in params:
            self.set_param(key, params[key])
