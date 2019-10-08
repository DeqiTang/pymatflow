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
        self.kptopt = 1
        self.ngkpt = [1, 1, 1]

        self.nshiftk = 1 
        #self.shiftk = np.zeros([self.nshiftk, 3])
        self.shiftk = np.full([self.nshiftk, 3], 0.5)
    
    def to_in(self, fout):
        # fout: a file stream for writing
        fout.write("# kpoints setting\n")
        #
        if self.kptopt == 1:
            fout.write("kptopt 1\n")
            fout.write("ngkpt %d %d %d\n" %(self.ngkpt[0], self.ngkpt[1], self.ngkpt[2]))
            fout.write("nshiftk %d\n" % self.nshiftk)
            fout.write("shiftk\n")
            for i in range(self.nshiftk):
                fout.write("%f %f %f\n" % (self.shiftk[i][0], self.shiftk[i][1], self.shiftk[i][2]))
        #
        if self.kptopt == 0:
            fout.write("kptopt 0\n")
            fout.write("nkpt \n")
            fout.write("kpt \n")
            fout.write("kptnrm \n")
            fout.write("wtk \n")
        #
        if self.kptopt == 2:
            pass
        #
        fout.write("\n")

    def set_kpoints(self, kptopt=1, ngkpt=[1, 1, 1]):
        self.kptopt = koptopt
        self.ngkpt = ngkpt

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
        fout.write("# electronic structure setting\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item, str(self.params[item])))
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
            if i in self.params.keys():
                if self.params[i] is not 0:
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
