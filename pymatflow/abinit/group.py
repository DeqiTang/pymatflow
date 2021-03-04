from pymatflow.variable import Variable


class AbinitVariableGroup:

    def __init__(self):
        self.params = {}
        self.n = 0
        self.status = True

    def to_string(self):
        self.sync_n()
        if False == self.status:
            return ""
        out = ""

        for item in self.params:
            if None == self.params[item].as_val():
                continue
            out += self.params[item].to_string(layout="same-line", indent="") + "\n"

        return out

    def set_param(self, key, value):
        self.remove(key)
        self.params[key] = Variable(key, value)
        self.params[key].n = self.n

    def contains(self, key):
        if key in self.params:
            return True
        else:
            return False

    def set_status(self, key, status):
        if False == self.contains(key):
            return
        else:
            self.params[key].status = status

    def remove(self, key):
        if key in self.params:
            del self.params[key]

    def clear(self):
        self.params.clear()

    def get(self, key, t=str, dim=0):
        """
        Note:
            return a 2D array of t type if dim == 2
            return a 1D array of t type if dim == 1
            return a scalar of t type if dim == 0
        """
        return self.params[key].as_val(t=t, dim=dim)

    def sync_n(self):
        """
        Note:
            synchronize Variable().n to self.n
        """
        for key in self.params:
            self.params[key].n = self.n

    def set_n(self, n):
        self.n = n
        self.sync_n()