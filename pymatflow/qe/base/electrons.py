"""
in control of &electrons /
"""
import sys

from pymatflow.qe.group import QeVariableGroup
from . import qe_variable_to_string
"""
usage:
"""

class QeElectrons(QeVariableGroup):
    """

    """
    def __init__(self):
        super().__init__()
        self.set_params({
                "electron_maxstep": None,
                "scf_must_converge": None,
                "conv_thr": None,
                "adaptive_thr": None,
                "conv_thr_init": None,
                "conv_thr_multi": None,
                "mixing_mode": None,
                "mixing_beta": None,
                "mixing_ndim": None,
                "mixing_fixed_ns": None,
                "diagonalization": None,
                "ortho_para": None,
                "diago_thr_init": None,
                "diago_cg_maxiter": None,
                "diago_david_ndim": None,
                "diago_full_acc": None,
                "efield": None,
                "efield_cart": None,
                "efield_phase": None,
                "startingpot": None,
                "startingwfc": None,
                "tqr": None,
                "real_space": None,
        })
    def to_in(self, fout):
        """
        ;param fout: a file stream for writing
        """
        fout.write(self.to_string())

    def to_string(self):
        out = ""        
        out += "&electrons\n"
        for item in self.params:
            if self.params[item].as_val() == None:
                continue            
            out += qe_variable_to_string(self.params[item])
            out += "\n"
        out += "/\n"
        out +=  "\n"
        return out

    def basic_setting(self):
        self.set_param("conv_thr", 1.0E-6)
        self.set_param("mixing_mode", "plain") # namely charge density Broyden mixing
        self.set_param("mixing_beta", 0.7E0) # mixing factor for self-consistency
        # number of iterations used in mixing scheme
        # if tight with memory, we can reduce it to 4
        self.set_param("mixing_ndim", 8)
        self.set_param("diagonalization", 'david')

