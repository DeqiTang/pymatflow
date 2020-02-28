
class abinit_misc:
    """
    """
    def __init__(self):
        self.params = {}
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
        for item in self.params:
            if self.params[item] is not None:
                input_str += "%s%s %s\n" % (item, n if n > 0 else "", str(self.params[item]))
                input_str += "\n"
        input_str += "\n\n"
        return input_str
        #


    def set_params(self, params):
        for item in params:
            self.params[item] = params[item]
