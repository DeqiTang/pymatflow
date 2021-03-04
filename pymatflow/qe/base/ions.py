"""
in control of &ions /
"""
import sys

from pymatflow.qe.group import QeVariableGroup
from . import qe_variable_to_string
"""
usage:
"""

class QeIons(QeVariableGroup):
    """

    """
    def __init__(self):
        super().__init__()
        self.set_params({
                "ion_dynamics": None,
                "ion_positions": None,
                "pot_extrapolation": None,
                "wfc_extrapolation": None,
                "remove_rigid_rot": None,
                "ion_temperature": None,
                "tempw": None,
                "tolp": None,
                "delta_t": None,
                "nraise": None,
                "refold_pos": None,
                "upscale": None,
                "bfgs_ndim": None,
                "trust_radius_max": None,
                "trust_radius_min": None,
                "trust_radius_ini": None,
                "w_1": None,
                "w_2": None,
        })

    def to_in(self, fout):
        """
        :param fout: a file stream for writing
        """
        fout.write(self.to_string())
    
    def to_string(self):
        out = ""
        out += "&ions\n"
        for item in self.params:
            if self.params[item].as_val() == None:
                continue            
            out += qe_variable_to_string(self.params[item])
            out += "\n"
        out += "/\n"
        out += "\n"
        return out

    def basic_setting(self, calc='relax'):
        """
            for different kind of running set different parameters
        """
        if calc == 'relax':
            self.set_param("ion_dynamics", "bfgs")
            #self.set_param("ion_temperature", 'not-controlled')
            #self.set_param("tempw", 300)
        elif calc == 'md':
            self.set_param("ion_dynamics", 'verlet')
            self.set_param("ion_temperature", 'not_controlled')
            self.set_param("tempw", 300)
        elif calc == 'vc-relax':
            self.set_param("ion_dynamics", 'bfgs')
            #self.set_param("ion_temperature", 'notcontrolled')
            #self.set_param("tempw", 300)
