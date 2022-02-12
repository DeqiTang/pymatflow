from pymatflow.abinit.group import AbinitVariableGroup

class AbinitMisc(AbinitVariableGroup):
    """
    """
    def __init__(self):
        super().__init__()
        #self.incharge = []
        self.status = True

    def to_string(self, n=0):
        """
        :return input_str is the string of all the set params
        """
        input_str = ""
        input_str += "# ============================\n"
        input_str += "# miscellaneous parameters\n"
        input_str += "# ============================\n"
        input_str += "\n"
        self.set_n(n)
        input_str += super().to_string()
        return input_str
        #