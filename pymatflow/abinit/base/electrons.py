"""
in control of electrons step related parameters
"""
import sys
import numpy as np

class Kpoints:
    """
    default behavior:
        use koptopt = 1 and ngkpt = [1, 1, 1]
    Note:
        I found the setting of kpoints in abinit is versatile
        and I haven't grasped the command of most of them.
        So at present only use the simplest way namely using
        kptopt == 1.
    Reference:
        https://docs.abinit.org/variables/basic/#kptopt
        https://docs.abinit.org/variables/gstate/#kptbounds
    """
    def __init__(self):
        self.params = {}
        self.incharge = [
            "kptopt", "ngkpt", "nshiftk", "shiftk",
            ]
        self.basic_setting()

    def basic_setting(self):
        self.params["kptopt"] = 1
        self.params["ngkpt"] = [1, 1, 1]
        #self.params["nshiftk"] = 1
        #self.params["shiftk"] = np.full([self.params["nshiftk"], 3], 0.0)
        #self.params["shiftk"] = np.zeros([self.params["nshiftk"], 3])

    def to_string(self, n=0):
        """
        :return input_str is the string of all the set params
        """
        input_str = ""
        input_str += "# kpoints setting\n"
        #
        if self.params["kptopt"] == 1:
            input_str += "kptopt%s 1\n\n" % (n if n > 0 else "")
            input_str += "ngkpt%s %d %d %d\n\n" %(n if n > 0 else "", self.params["ngkpt"][0], self.params["ngkpt"][1], self.params["ngkpt"][2])
            #input_str += "nshiftk %d\n\n" % self.params["nshiftk"]
            #input_str += "shiftk\n"
            #for i in range(self.params["nshiftk"]):
            #    input_str += "%f %f %f\n" % (self.params["shiftk"][i][0], self.params["shiftk"][i][1], self.params["shiftk"][i][2])

            #input_str += "istwfk 1\n") # for rf
        #
        if self.params["kptopt"] == 0:
            input_str += "kptopt%s 0\n\n" % (n if n > 0 else "")
            input_str += "nkpt%s \n\n" % (n if n > 0 else "")
            input_str += "kpt%s \n\n" % (n if n > 0 else "")
            input_str += "kptnrm%s \n\n" % (n if n > 0 else "")
            input_str += "wtk%s \n\n" % (n if n > 0 else "")
            input_str += "istwfk%s 1\n" % (n if n > 0 else "")
        #
        if self.params["kptopt"] == 2:
            input_str += "kptopt%s 2\n" % (n if n > 0 else "")
            input_str += "ngkpt%s %d %d %d\n\n" %(n if n > 0 else "", self.params["ngkpt"][0], self.params["ngkpt"][1], self.params["ngkpt"][2])
            #input_str += "istwfk 1\n") # for rf
        if self.params["kptopt"] == 3:
            # typically for rf calculation
            input_str += "kptopt%s %d\n" % (n if n > 0 else "", self.params["kptopt"])
            input_str += "ngkpt%s %d %d %d\n\n" %(n if n > 0 else "", self.params["ngkpt"][0], self.params["ngkpt"][1], self.params["ngkpt"][2])
            #input_str += "nshiftk %d\n\n" % self.params["nshiftk"]
            #input_str += "shiftk\n"
            #for i in range(self.params["nshiftk"]):
            #    input_str += "%f %f %f\n" % (self.params["shiftk"][i][0], self.params["shiftk"][i][1], self.params["shiftk"][i][2])
        #
        if self.params["kptopt"] < 0:
            # for band structure calculation using kptbounds and ndivk(ndivsm), iscf must equal to -2
            # for format of self.params["kptbounds"] see self.set_band()
            input_str += "kptopt%s %d\n" % (n if n > 0 else "", self.params["kptopt"])
            input_str += "ndivk\n"
            for i in range(len(self.kpath) - 1):
                if self.kpath[i][4] != "|":
                    input_str += "%d " % (self.kpath[i][4])
                else:
                    input_str += "%d " % (1)
            input_str += "\n"
            input_str += "kptbounds%s\n" % (n if n > 0 else "")
            for i in range(len(self.kpath)):
                input_str += "%f %f %f #%s %s\n" % (
                    self.kpath[i][0],
                    self.kpath[i][1],
                    self.kpath[i][2],
                    self.kpath[i][3],
                    self.kpath[i][4],
                    )
        # end
        #
        input_str += "\n"

        return input_str


    def set_band(self, kpath):
        """
        :parma kpath:
                the high symmetry k point path used in bands structure calculation
                in format like this:

                [[kx, ky, kz, label, connect_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', 15], ...]

                if connect_indicator in a kpoint is an integer, then it will connect to the following point
                through the number of kpoints defined by connect_indicator.

                if connect_indicator in a kpoint is '|', then it will not connect to the following point,
        """
        self.kpath = kpath
        self.params["kptopt"] = - (len(self.kpath) - 1)
        #

    def set_params(self, kpoints):
        for item in kpoints:
            self.params[item] = kpoints[item]


class AbinitElectrons:
    """
    """
    def __init__(self):
        self.params = {}
        self.incharge = [
            "ecut", "ixc", "nstep", "diemac", "iscf",
            'toldfe', 'tolwfr', 'toldff', 'tolrff', 'tolvrs',
            "occopt", "nband", "occ", "wtk",
            "prtden", "prtdos",
            ]
        self.status = True
        self.kpoints = Kpoints()

    def to_string(self, n=0):
        """
        :return input_str is the string of all the set params
        """
        input_str = ""
        # ------------
        # 检查输入参数
        self.check_all_params()
        # ---------------
        # 检查输入参数结束
        input_str += "# =====================================\n"
        input_str += "# electronic structure related setting\n"
        input_str += "# =====================================\n"
        input_str += "\n"
        for item in self.params:
            if self.params[item] is not None:
                input_str +="%s%s %s\n" % (item, n if n > 0 else "", str(self.params[item]))
                input_str += "\n"
        #fout.write("\n")
        # 输入k点
        input_str += self.kpoints.to_string(n=n)
        #
        input_str += "\n"

        return input_str

    def check_all_params(self):
        # this function is responsible for calling other
        # check function to control the overall parameters
        # checking.
        self.check_scf_criteria()
        self.check_kpoints()


    def use_tol(self, tol, value):
        """
        :param tol: one of ['toldfe', 'tolwfr', 'toldff', 'tolrff', 'tolvrs'].
        :param value: set value for the choosen tol criteria

        Reference: see https://docs.abinit.org/topics/SCFControl/
        """
        tols = ['toldfe', 'tolwfr', 'toldff', 'tolrff', 'tolvrs']
        for item in tols:
            if item != tol:
                self.params[item] = None
        self.params[tol] = value



    def check_scf_criteria(self):
        """
            there are six kinds of criteria to judge convergence of
            scf in abinit, but one and only one of them can be used
            in a run. here we make a check
        """
        tols = ['toldfe', 'tolwfr', 'toldff', 'tolrff', 'tolvrs']
        nonzeros = 0
        for i in tols:
            if i in self.params.keys() and self.params[i] is not None:
                if self.params[i] != 0.0:
                    nonzeros += 1
        if nonzeros == 1:
            return True
        else:
            print("========================================\n")
            print("                WARNING !!!\n")
            print("========================================\n")
            print("you must set one and only one of variables\n")
            print("below to differ from zero.\n")
            print("[toldfe, tolwfr, toldff, tolrff, tolvrs]\n")
            #print(nonzeros)
            sys.exit(1)

    def check_kpoints(self):
        # there should be no kpoints related setting in self.params
        tmp = [
                'kptopt', 'ngkpt', 'istwfk', 'kpt', 'kptbounds',
                'kptnrm', 'kptrlatt', 'kptrlen', 'ndivk', 'ndivsm',
                'nkpath', 'nkpt', 'nshiftk', 'prtkpt', 'shiftk', 'wtk'
              ]
        for i in tmp:
            if i in self.params:
                print("========================================\n")
                print("         Warning !!!\n")
                print("========================================\n")
                print("do not set kpoints through abinit_electrons.params\n")
                print("as you should do it via abinit_electrons.kpoints.\n")
                sys.exit(1)
        #

    def dft_plus_u(self):
        # 需要使用paw赝势, 在files文件中进行设置
        self.params["usepawu"] = 1
        self.params["pawecutdg"] = 30
        self.params["lpawu"] = [-1, -1, -1, 2]
        self.params["upawu"] = [4.0, 4.0, 4.0, 4.0]
        self.params["jpawu"] = [0.8, 0.8, 0.8, 0.8]

    def basic_setting(self):
        self.params["ecut"] = 15 #15
        #self.params["pawecutdg"] = 50
        self.params["occopt"] = 3  # fermi dirac smearing of occupation
        self.params["nstep"] = 100
        self.params["diemac"] = 2.0
        self.use_tol(tol="tolvrs", value = 1.0e-18)
        self.params["ixc"] = 11

    def set_params(self, params):
        for item in params:
            self.params[item] = params[item]

    def set_scf_nscf(self, mode="scf"):
        if mode == "scf":
            if "usepaw" in self.params and self.params["usepaw"] == 1:
                self.params["iscf"] = 17
            else:
                self.params["iscf"] = 7
            self.params["prtden"] = 1
        if mode == "nscf":
            # -3 is good for band structure calculation
            # -2 is good for dos calculation
            self.params["iscf"] = -3
            self.params["nstep"] = 0
