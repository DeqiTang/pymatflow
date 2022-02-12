"""
in control of &cell /
"""
import sys

from pymatflow.qe.group import QeVariableGroup
from . import qe_variable_to_string

"""
usage:
"""

class QeCell(QeVariableGroup):
    """

    """
    def __init__(self):
        super().__init__()

    def to_in(self, fout):
        fout.write(self.to_string())

    def to_string(self):
        # fout: a file stream for writing
        out = ""
        out += "&cell\n"
        for item in self.params:
            if self.params[item].as_val() == None:
                continue            
            out += qe_variable_to_string(self.params[item])
            out += "\n"
        out += "/\n"
        out += "\n"

