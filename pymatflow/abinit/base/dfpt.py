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
        self.status = False
        #self.kpoints = kpoints()

    def to_string(self, n=0):
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
                    input_str += "rfdir%s %d %d %d\n" % (n if n > 0 else "", self.params[item][0], self.params[item][1], self.params[item][2])
                elif item == "rfatpol":
                    input_str += "rfatpol%s %d %d\n" % (n if n > 0 else "", self.params["rfatpol"][0], self.params["rfatpol"][1])
                    pass
                elif item == "qpt":
                    input_str += "qpt%s %f %f %f\n" % (n if n > 0 else "", self.params[item][0], self.params[item][1], self.params[item][2])
                elif item == "ngqpt":
                    input_str += "ngqpt%s %d %d %d\n" % (n if n > 0 else "", self.params[item][0], self.params[item][1], self.params[item][2])
                else:
                    input_str += "%s%s %s\n" % (item, n if n > 0 else "", str(self.params[item]))
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
