"""

"""

class hsd_block:
    """
    a hsd_block may have one of the two types: property | method
    """
    def __init__(self, name, block_type, val=None, level=0):
        self.name = name
        self.level = level
        self.block_type = block_type
        self.val = val
        self.status = False

        self.scalar = {}
        self.list_of_scalar = {}
        self.method = {}
        self.list_of_property = {}

    def to_string(self):
        
        if self.status == False:
            return ""

        out = ""
        indent = " " * self.level
        
        if self.block_type == "method":
            out += indent + "%s = %s {\n" % (self.name, self.val if self.val != None else "")
        else:
            out += indent + "%s = {\n" %(self.name)

        for item in self.scalar:
            out += indent + "%s = %s\n" % (item, self.scalar[item])
        for item in self.list_of_scalar:
            if item == "":
                out += indent
                for val in self.list_of_scalar[item]:
                    out += " %s" % val
                out += "\n"
            else:
                out += indent + "%s ="
                for val in self.list_of_scalar[item]:
                    out += " %s" % val
                out += "\n"
        for item in self.method:
            out += self.method[item].to_string()
        for item in self.list_of_property:
            out += self.list_of_property[item].to_string()
        
        out += indent + "}\n"

        return out
        
