
class abinit_misc:
    """
    """
    def __init__(self):
        self.params = {}
        #self.incharge = []

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("# ============================\n")
        fout.write("# miscellaneous parameters\n")
        fout.write("# ============================\n")
        fout.write("\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item ,str(self.params[item])))
                fout.write("\n")
        fout.write("\n\n")
        #

    def set_params(self, params):
        for item in params:
            self.params[item] = params[item]
