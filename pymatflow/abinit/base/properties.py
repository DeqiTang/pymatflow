"""
in control of properties calculation  related parameters
"""
from pymatflow.abinit.group import AbinitVariableGroup


class AbinitProperties(AbinitVariableGroup):
    """
    """
    def __init__(self):
        super().__init__()
        self.incharge = []
        self.status = True

    def to_string(self, n=0):
        """
        :return input_str is the string of all the set params
        """
        input_str = ""
        input_str += "# ======================================\n"
        input_str += "# properties calculation related setting\n"
        input_str += "# ======================================\n"
        input_str += "\n"
        self.set_n(n)
        input_str += super().to_string()
        input_str += "\n"
        input_str += "\n"
        return input_str
        #


    def get_option(self, option=[]):
        """
        :param  1: dos
                3: elf
                4:
        """
        if len(option) == 0:
            return
        if 1 in option:
            self.do_dos()
        if 2 in option:
            self.do_bands()
        if 3 in option:
            self.do_elf()
        if 4 in option:
            pass

    def do_dos(self):
        self.set_param("prtdos", 2)

    def do_elf(self):
        self.set_param("prtelf", 1)

    def berry_phase(self):
        self.set_param("berryopt", -1)
        self.set_param("nberry", 8)
        self.set_params("dberry", [1, 1, 1, 1])
        self.set_param("rfdir", [1, 1, 1])

    def piezoelectric(self):
        self.set_param("piezoflag", 3)
