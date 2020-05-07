

class execution:
    """
    """
    def __init__(self):
        self.params = {}

    def to_string(self):
        out  = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out