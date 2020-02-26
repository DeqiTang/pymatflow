"""
in control of DFPT related parameters
"""
import sys
import numpy as np


class abinit_dfpt:
    """
    Reference:
        https://docs.abinit.org/guide/respfn/
        https://docs.abinit.org/tutorial/rf1/
        https://docs.abinit.org/tutorial/rf2/
        https://docs.abinit.org/tutorial/elastic/
        https://docs.abinit.org/tutorial/ffield/
        https://docs.abinit.org/tutorial/nlo/
        https://docs.abinit.org/topics/DFPT/
        https://docs.abinit.org/variables/dfpt/
    """
    def __init__(self):
        self.params = {}
        self.incharge = [
            "rfatpol", "rfddk", "rfelfd", "rfphon", "rfstrs", "rfdir", "rfmagn",
            "rfmeth", "asr", "dfpt_sciss",
            ]
        #self.kpoints = kpoints()

    def to_input(self, fout):
        # fout: a file stream for writing
        # ------------
        # 检查输入参数
        #self.check_all_params()
        # ---------------
        # 检查输入参数结束
        fout.write("# =====================================\n")
        fout.write("# DFPT related setting\n")
        fout.write("# =====================================\n")
        fout.write("\n")
        for item in self.params:
            if self.params[item] is not None:
                if item == "rfdir":
                    fout.write("rfdir %d %d %d\n" % (self.params[item][0], self.params[item][1], self.params[item][2]))
                elif item == "rfatpol":
                    fout.write("rfatpol %d %d\n" % (self.params["rfatpol"][0], self.params["rfatpol"][1]))
                    pass
                elif item == "qpt":
                    fout.write("qpt %f %f %f\n" % (self.params[item][0], self.params[item][1], self.params[item][2]))
                elif item == "ngqpt":
                    fout.write("ngqpt %d %d %d\n" % (self.params[item][0], self.params[item][1], self.params[item][2]))
                else:
                    fout.write("%s %s\n" % (item, str(self.params[item])))
                fout.write("\n")
        fout.write("\n")

    def to_string(self):
        """
        :return input_str is the string of all the set params
        """
        # ------------
        # 检查输入参数
        #self.check_all_params()
        # ---------------
        # 检查输入参数结束
        input_str = ""
        input_str += "# =====================================\n"
        input_str += "# DFPT related setting\n"
        input_str += "# =====================================\n"
        input_str += "\n"
        for item in self.params:
            if self.params[item] is not None:
                if item == "rfdir":
                    input_str += "rfdir %d %d %d\n" % (self.params[item][0], self.params[item][1], self.params[item][2])
                elif item == "rfatpol":
                    input_str += "rfatpol %d %d\n" % (self.params["rfatpol"][0], self.params["rfatpol"][1])
                    pass
                elif item == "qpt":
                    input_str += "qpt %f %f %f\n" % (self.params[item][0], self.params[item][1], self.params[item][2])
                elif item == "ngqpt":
                    input_str += "ngqpt %d %d %d\n" % (self.params[item][0], self.params[item][1], self.params[item][2])
                else:
                    input_str += "%s %s\n" % (item, str(self.params[item]))
                input_str += "\n"
        input_str += "\n"

        return input_str


    def basic_setting(self):
        self.params["rfdir"] = [1, 1, 1]
        #self.params["nqpt"] = 1
        #self.params["qpt"] = [0, 0, 0]

    def set_params(self, params):
        for item in params:
            self.params[item] = params[item]
