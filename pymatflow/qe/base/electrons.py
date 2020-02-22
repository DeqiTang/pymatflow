"""
in control of &electrons /
"""
import sys


"""
usage:
"""

class qe_electrons:
    """

    """
    def __init__(self):
        self.params = {
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
                }
    def to_in(self, fout):
        """
        ;param fout: a file stream for writing
        """
        fout.write("&electrons\n")
        for item in self.params:
            if self.params[item] is not None:
                if type(self.params[item]) is str:
                    if self.params[item] == ".true." or self.params[item] == ".false.":
                        fout.write("%s = %s\n" % (item, str(self.params[item])))
                    else:
                        fout.write("%s = '%s'\n" % (item, str(self.params[item])))
                else:
                    fout.write("%s = %s\n" % (item, str(self.params[item])))
        fout.write("/\n")
        fout.write("\n")

    def basic_setting(self):
        self.params["conv_thr"] = 1.0E-6
        self.params["mixing_mode"] = "plain" # namely charge density Broyden mixing
        self.params["mixing_beta"] = 0.7E0 # mixing factor for self-consistency
        # number of iterations used in mixing scheme
        # if tight with memory, we can reduce it to 4
        self.params["mixing_ndim"] = 8
        self.params["diagonalization"] = 'david'

    def set_params(self, params):
        """
        :param params: a dict storing the parameters and values
        """
        for item in params:
            self.params[item] = params[item]
