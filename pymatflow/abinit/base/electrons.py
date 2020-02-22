"""
in control of electrons step related parameters
"""
import sys
import numpy as np

class kpoints:
    """
    default behavior:
        use koptopt == 1 and ngkpt = [1, 1, 1]
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
        self.params["nshiftk"] = 1
        self.params["shiftk"] = np.full([self.params["nshiftk"], 3], 0.5)
        #self.params["shiftk"] = np.zeros([self.params["nshiftk"], 3])


    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("# kpoints setting\n")
        #
        if self.params["kptopt"] == 1:
            fout.write("kptopt 1\n\n")
            fout.write("ngkpt %d %d %d\n\n" %(self.params["ngkpt"][0], self.params["ngkpt"][1], self.params["ngkpt"][2]))
            fout.write("nshiftk %d\n\n" % self.params["nshiftk"])
            fout.write("shiftk\n")
            for i in range(self.params["nshiftk"]):
                fout.write("%f %f %f\n" % (self.params["shiftk"][i][0], self.params["shiftk"][i][1], self.params["shiftk"][i][2]))

            #fout.write("istwfk 1\n") # for rf
        #
        if self.params["kptopt"] == 0:
            fout.write("kptopt 0\n\n")
            fout.write("nkpt \n\n")
            fout.write("kpt \n\n")
            fout.write("kptnrm \n\n")
            fout.write("wtk \n\n")
            fout.write("istwfk 1\n") # for rf
        #
        if self.params["kptopt"] == 2:
            fout.write("kptopt 2\n")
            fout.write("ngkpt %d %d %d\n\n" %(self.params["ngkpt"][0], self.params["ngkpt"][1], self.params["ngkpt"][2]))
            #fout.write("istwfk 1\n") # for rf
        if self.params["kptopt"] == 3:
            # typically for rf calculation
            fout.write("kptopt %d\n" % self.params["kptopt"])
            fout.write("ngkpt %d %d %d\n\n" %(self.params["ngkpt"][0], self.params["ngkpt"][1], self.params["ngkpt"][2]))
            fout.write("nshiftk %d\n\n" % self.params["nshiftk"])
            fout.write("shiftk\n")
            for i in range(self.params["nshiftk"]):
                fout.write("%f %f %f\n" % (self.params["shiftk"][i][0], self.params["shiftk"][i][1], self.params["shiftk"][i][2]))
        #
        if self.params["kptopt"] < 0:
            # for band structure calculation using kptbounds and ndivk(ndivsm), iscf must equal to -2
            # for format of self.params["kptbounds"] see self.set_band()
            fout.write("kptopt %d\n" % self.params["kptopt"])
            fout.write("ndivsm %d\n" % 10)
            fout.write("kptbounds\n")
            for i in range(len(self.params["kptbounds"])):
                fout.write("%f %f %f #%s\n" % (
                    self.params["kptbounds"][i][0],
                    self.params["kptbounds"][i][1],
                    self.params["kptbounds"][i][2],
                    self.params["kptbounds"][i][3],
                    ))
        # end
        #
        fout.write("\n")

    def set_band(self, kptbounds=None):
        """
            self.params["kptbounds"]:
                the high symmetry k point path used in bands structure calculation
                in format like this:

                [[kx, ky, kz, label, connect_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', 15], ...]

                if connect_indicator in a kpoint is an integer, then it will connect to the following point
                through the number of kpoints defined by connect_indicator.

                if connect_indicator in a kpoint is '|', then it will not connect to the following point,
        """
        self.params["kptbounds"] = kptbounds
        self.params["kptopt"] = - (len(self.params["kptbounds"]) - 1)
        #

    def set_params(self, kpoints):
        for item in kpoints:
            self.params[item] = kpoints[item]


class abinit_electrons:
    """
    """
    def __init__(self):
        self.params = {}
        self.incharge = [
            "ecut", "ixc", "nstep", "toldfe", "diemac"
            ]
        self.kpoints = kpoints()

    def to_input(self, fout):
        # fout: a file stream for writing
        # ------------
        # 检查输入参数
        self.check_all_params()
        # ---------------
        # 检查输入参数结束
        fout.write("# =====================================\n")
        fout.write("# electronic structure related setting\n")
        fout.write("# =====================================\n")
        fout.write("\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item, str(self.params[item])))
                fout.write("\n")
        #fout.write("\n")
        # 输入k点
        self.kpoints.to_input(fout)
        #
        fout.write("\n")


    def check_all_params(self):
        # this function is responsible for calling other
        # check function to control the overall parameters
        # checking.
        self.check_scf_criteria()
        self.check_kpoints()

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
        self.params["lpawu"] = '-1 -1 -1 2'
        self.params["upawu"] = "4.0 4.0 4.0 4.0"
        self.params["jpawu"] = "0.8 0.8 0.8 0.8"

    def basic_setting(self):
        self.params["ecut"] = 15 #15
        #self.params["pawecutdg"] = 50
        self.params["occopt"] = 3  # fermi dirac smearing of occupation
        self.params["nstep"] = 100
        self.params["diemac"] = 2.0
        #self.params["toldfe"] = 1.0e-6
        self.params["tolvrs"] = 1.0e-18
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
