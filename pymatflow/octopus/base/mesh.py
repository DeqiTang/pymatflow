

class mesh:
    """
    """
    def __init__(self):
        self.params = {}
                
        self.kpoints_option = "mp"
        self.kpoints_mp = [1, 1, 1, 0, 0, 0]

    def to_string(self):
        out  = ""
        for item in self.params:
            out += "%s = %s\n" % (item, self.params[item])
            out += "\n"
        return out

    def write_kpoints(self, fout):
        """
        :param fout: a file stream for writing
        """
        if self.kpoints_option == "mp":
            fout.write("\%KPointsGrid\n")
            fout.write("%d | %d | %d\n" % (
                self.kpoints_mp[0],
                self.kpoints_mp[1],
                self.kpoints_mp[2],
                ))
            fout.write("\%\n")                
        elif self.kpoints_option == "kpath":
            # there is a trick:
            # when self.kpath[i][4] == "|"
            # we set the number of k point to connect to the next high symmetry kpoint to 0
            # this is very fantastic !!!
            fout.write("\%KPointsPath\n")
            for i in range(len(self.kpath)-2):
                if self.kpath[i][4] == "|":
                    fout.write("0 | ")
                else:
                    fout.write("%d | " % self.kpath[i][4])
                fout.write("%d\n" % self.kpath[i][-2])
            #
            for i in range(len(self.kpath)):
                fout.write("%f | %f | %f #%s\n" % sefl.kpath[i][3])
            #
            #fout.write("KPointsUseSymmetries = no\n")
            fout.write("\%\n")      

    def set_kpoints(self, kpoints_mp=[1, 1, 1, 0, 0, 0], option="mp", kpath=None):
        """
        :param kpath: the high symmetry k point path used in bands structure calculation
            in format like this:

            [[kx, ky, kz, label, connect_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', 15], ...]

            if connect_indicator in a kpoint is an integer, then it will connect to the following point
            through the number of kpoints defined by connect_indicator.

            if connect_indicator in a kpoint is '|', then it will not connect to the following point,
        TODO:
        Note:
            "mp" means Monkhorst-Pack scheme
        """
        if option == "mp":
            self.kpoints_option = option
            self.kpoints_mp = kpoints_mp
            return
        if option == "kpath":
            self.kpoints_option = option
            self.kpath = kpath
            return