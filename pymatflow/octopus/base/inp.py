
from pymatflow.octopus.base.calculation_modes import calculation_modes
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
        self.calculation_modes = calculation_modes()
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
        out += self.calculation_mode.to_string()
        out += self.execution.to_string()
        out += self.hamiltonian.to_string()
        out += self.linear_response.to_string()
        out += self.math.to_string()
        out += self.mesh.to_string()
        out += self.output.to_string()
        out += self.scf.to_string()
        out += self.states.to_string()
        out += self.system.to_string()
        out += self.time_dependent.to_string()
        out += self.utilities.to_string()

    def set_runtype(self, runtype="static"):
        """
        self.runtype: static | opt | md | dfpt | phonon | neb
        """
        self.runtype = runtype
        self.set_calculation_mode_default(mode="gs")
        
    def set_calculation_mode_default(self, mode="gs"):
        """
        self.calculation_mode: gs | unocc | td | go | opt_control | em_resp | casida | vdw | vib_modes | one_shot | kdotp | dummy | invert_ks | recipe
        Note:
            set the default parameters for the chosen CalculationMode
        """
        self.calculation_mode.set_default(mode=mode)

    def set_params(self, params):
        """
        """
        for item in params:
            if item.split("/")[0] == "Calculation Modes":
                self.calculation_modes.set_params({item: params[item]})
            elif item.split("/")[0] == "Execution":
                self.execution.set_params({item: params[item]})
            elif item.split("/")[0] == "Hamiltonian":
                self.hamiltonian.set_params({item: params[item]})
            elif item.split("/")[0] == "LinearResponse":
                self.linear_response.set_params({item: params[item]})
            elif item.split("/")[0] == "Math":
                self.math.set_params({item: params[item]})
            elif item.split("/")[0] == "Mesh":
                self.mesh.set_params({item: params[item]})
            elif item.split("/")[0] == "Output":
                self.output.set_params({item: params[item]})
            elif item.split("/")[0] == "SCF":
                self.scf.set_params({item: params[item]})
            elif item.split("/")[0] == "States":
                self.states.set_params({item: params[item]})
            elif item.split("/")[0] == "System":
                self.system.set_params({item: params[item]})
            elif item.split("/")[0] == "Time-Dependent":
                self.time_dependent.set_params({item: params[item]})
            elif item.split("/")[0] == "Utilitites":
                self.utilities.set_params({item: params[item]})
            else:
                pass