
from pymatflow.octopus.base.calculation_mode import calculation_mode
from pymatflow.octopus.base.execution import execution
from pymatflow.octopus.base.hamiltonian import hamiltonian
from pymatflow.octopus.base.linear_response import linear_response
from pymatflow.octopus.base.math import math
from pymatflow.octopus.base.mesh import mesh
from pymatflow.octopus.base.output import output
from pymatflow.octopus.base.scf import scf
from pymatflow.octopus.base.states import states
from pymatflow.octopus.base.system import system 
from pymatflow.octopus.base.time_dependent import time_dependent
from pymatflow.octopus.base.utilities import utilities


class inp:
    """
    """
    def __init__(self):
        self.calculation_mode = calculation_mode()
        self.execution = execution()
        self.hamiltonian = hamiltonian()
        self.linear_response = linear_response()
        self.math = math()
        self.mesh = mesh()
        self.output = output()
        self.scf = scf()
        self.states = states()
        self.system = system()
        self.time_dependent = time_dependent()
        self.utilities = utilities()

    def to_inp(self, fout):
        fout.write(self.calculation_mode.to_string())
        fout.write(self.execution.to_string())
        fout.write(self.hamiltonian.to_string())
        fout.write(self.linear_response.to_string())
        fout.write(self.math.to_string())
        fout.write(self.mesh.to_string())
        fout.write(self.output.to_string())
        fout.write(self.scf.to_string())
        fout.write(self.states.to_string())
        fout.write(self.system.to_string())
        fout.write(self.time_dependent.to_string())
        fout.write(self.utilities.to_string())

    def to_string(self):
        out = ""
        out += self.calculation_mode.to_string())
        out += self.execution.to_string())
        out += self.hamiltonian.to_string())
        out += self.linear_response.to_string())
        out += self.math.to_string())
        out += self.mesh.to_string())
        out += self.output.to_string())
        out += self.scf.to_string())
        out += self.states.to_string())
        out += self.system.to_string())
        out += self.time_dependent.to_string())
        out += self.utilities.to_string())
