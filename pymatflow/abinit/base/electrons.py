#!/usr/bin/env python
# _*_ coding: utf-8 _*_

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
    """
    def __init__(self):
        self.params = {
                "kptopt": None,
                "ngkpt": None,
                "nshiftk": None,
                "shiftk": None,
                }
        self.basic_setting()

    def basic_setting(self):
        self.params["kptopt"] = 1
        self.params["ngkpt"] = [1, 1, 1]
        self.params["nshiftk"] = 1
        self.params["shiftk"] = np.full([self.params["nshiftk"], 3], 0.5)
        #self.params["shiftk"] = np.zeros([self.params["nshiftk"], 3])

    
    def to_in(self, fout):
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
            fout.write("kptopt %d\n" % self.params["kptopt"])
            fout.write("ndivsm %d\n" % 10)
            fout.write("kptbounds\n")
            point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][0][0]]
            fout.write("%f %f %f #%s\n" % (point[0], point[1], point[2], self.kpoints_seekpath["path"][0][0]))
            point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][0][1]]
            fout.write("%f %f %f #%s\n" % (point[0], point[1], point[2], self.kpoints_seekpath["path"][0][1]))
            for i in range(1, len(self.kpoints_seekpath["path"])):
                if self.kpoints_seekpath["path"][i][0] == self.kpoints_seekpath["path"][i-1][1]:
                    point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][i][1]]
                    fout.write("%f %f %f #%s\n" % (point[0], point[1], point[2], self.kpoints_seekpath["path"][i][1]))
                else:
                    point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][i][0]]
                    fout.write("%f %f %f #%s\n" % (point[0], point[1], point[2], self.kpoints_seekpath["path"][i][0]))
                    point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][i][1]]
                    fout.write("%f %f %f #%s\n" % (point[0], point[1], point[2], self.kpoints_seekpath["path"][i][1]))
        # end
        #
        fout.write("\n")

    def set_band(self, system):
        """
        system is an instance of abinit.base.system.abinit_system
        """
        import seekpath
        lattice = [system.xyz.cell[0:3], system.xyz.cell[3:6], system.xyz.cell[6:9]]
        positions = []
        numbers = []
        a = np.sqrt(system.xyz.cell[0]**2 + system.xyz.cell[1]**2 + system.xyz.cell[2]**2)
        b = np.sqrt(system.xyz.cell[3]**2 + system.xyz.cell[4]**2 + system.xyz.cell[5]**2)
        c = np.sqrt(system.xyz.cell[6]**2 + system.xyz.cell[7]**2 + system.xyz.cell[8]**2)
        for atom in system.xyz.atoms:
            positions.append([atom.x / a, atom.y / b, atom.z / c])
            numbers.append(system.xyz.specie_labels[atom.name])
        structure = (lattice, positions, numbers)
        self.kpoints_seekpath = seekpath.get_path(structure)
        nks = 2
        for i in range(1, len(self.kpoints_seekpath["path"])):
            if self.kpoints_seekpath["path"][i][0] == self.kpoints_seekpath["path"][i-1][1]:
                nks = nks + 1
            else:
                nks = nks + 2
        self.params["kptopt"] = - (nks - 1)
        #

    def set_params(self, kpoints):
        for item in kpoints:
            self.params[item] = kpoints[item]

class abinit_electrons:
    """
    """
    def __init__(self):
        self.params = {
                "ecut": None,
                "ixc": None,
                "nstep": None,
                "toldfe": None,
                "diemac": None,
                }
        self.kpoints = kpoints()
                
    def to_in(self, fout):
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
        self.kpoints.to_in(fout)
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
        abinit中scf收敛盘踞有六种形式, 但是一次只能使用一种.
        这里对其进行检查
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


