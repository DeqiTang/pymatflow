from pymatflow.variable import Variable


class VariableGroup:
    """
    Note:
        derivative clas of VariableGroup should implement to_string() function themselves.
    """
    def __init__(self):
        self.params = {}
        self.status = True

    def set_param(self, key, value, unit=None):
        if key in self.params:
            self.params[key].set(key=key, value=value, unit=unit)
        else:
            self.params[key] = Variable(key, value, unit)

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
