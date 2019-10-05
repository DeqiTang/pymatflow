#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from emuhelper.abinit.base.electrons import abinit_electrons
from emuhelper.abinit.base.ions import abinit_ions
from emuhelper.abinit.base.system import abinit_system

"""
Usage:
Note:
    参考: https://atztogo.github.io/phonopy/abinit.html
"""

        
class phonopy:
    """
    """
    def __init__(self, xyz_f):
        self.system = abinit_system(xyz_f)
        self.electrons = abinit_electrons()
        self.ions = abinit_ions()

        self.electrons.params["ecut"] = 50
        self.electrons.params["kptopt"] = 1
        self.electrons.params["ngkpt"] = "1 1 1"
        self.electrons.params["occopt"] = 3  # fermi dirac smearing of occupation
        self.electrons.params["nstep"] = 100
        self.electrons.params["diemac"] = 2.0
        self.electrons.params["toldfe"] = 1.0e-6

        self.supercell_n = "1 1 1"
        self.head_inp_name = "head-phonon.in"
        self.pos_inp_name = "pos.in"

    def gen_input(self, directory="tmp-abinit-phonopy"):
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        os.system("cp *.psp8 %s/" % directory)

        cutoff = 40
        supercell_n = self.supercell_n

        
        os.chdir(directory)

        head_inp_name = self.head_inp_name
        with open(head_inp_name, 'w') as fout:
            self.electrons.to_in(fout)
            self.ions.to_in(fout)

        pos_inp_name = self.pos_inp_name
        with open(pos_inp_name, 'w') as fout:
            self.system.to_in(fout)

        os.system("phonopy --abinit -d --dim='%s' -c %s" % (supercell_n, pos_inp_name))
        
        disps = self.get_disps("./")
        
        for disp in disps:
            os.mkdir("disp-%s" % disp)
            os.system("cp ./*.psp8 ./disp-%s/" % (disp))
            # 解释一下下面进行的几个操作的意义:
            # 主要是我的脚本中XYZ在获取xyz文件中的结构信息时, 对元素装入self.specie_labels时进行了处理, 使得
            # 原子序数较小的原子对应的编号也更小. 这个编号顺序使得其构建abinit输入文件中的结构信息时对znucl
            # 后面的内容是从原子序数从小到大的顺序来设置的. 这样typat的设置会结合znucl中的顺序以及原子坐标顺序
            # 来给出指定. 这本来是没有问题的, 但是使用phonopy处理得到各个超胞时, 它没有保持我之前的znucl中按
            # 原子序数从小到大的顺序来设置, 但是原子坐标顺序没变, 不过typat对应也变了, 从结构的角度来说是正确的,
            # 不过这样一来znucl和typat变化使得xxx.files中给出赝势文件的顺序与其不一致.
            # 而我习惯于znucl原子序数从小到大, 因此下面的操作我就是想了一个办法, 将phonopy产生的超胞信息只保留
            # 原子坐标那一块儿, 然后使用之前为了创建超胞用的原始坐标文件中的znucl以及typat设置, 这样对应的结构
            # 也是没有错的, 而且原本xxx.files文件中给出元素赝势的顺序业余输入文件中znucl的一致, 便可以成功进行后续
            # 计算了.
            os.system("cat %s > ./disp-%s/supercell-%s.in" % (head_inp_name, disp, disp))
            os.system("head %s -n -%d >> ./disp-%s/supercell-%s.in" % (pos_inp_name, (self.system.xyz.natom + 2), disp, disp))
            os.system("tail supercell-%s.in -n %d >> ./disp-%s/supercell-%s.in" % (disp, (self.system.xyz.natom + 1), disp, disp))
            with open("./disp-%s/supercell-%s.files" % (disp, disp), 'w') as fout:
                fout.write("supercell-%s.in\n" % disp)
                fout.write("supercell-%s.out\n" % disp)
                fout.write("supercell-%si\n" % disp)
                fout.write("supercell-%so\n" % disp)
                fout.write("temp\n")
                for element in self.system.xyz.specie_labels:
                    fout.write("%s\n" % (element + ".psp8"))
        os.chdir("../")

    def run(self, directory="tmp-abinit-phonopy"):
        # run the simulation
        os.chdir(directory)
        disps = self.get_disps("./")
        for disp in disps:
            os.chdir("disp-%s" % disp)
            os.system("abinit < supercell-%s.files" % disp)
            os.chdir("../")
        os.chdir("../")
    
    def analysis(self, directory="tmp-abinit-phonopy"):
        # analyse the result
        os.chdir(directory)
        disps = self.get_disps("./")
        supercell_n = self.supercell_n
        pos_inp_name = self.pos_inp_name
        outs = ""
        for disp in disps:
            outs += "disp-%s/supercell-%s.out " % (disp, disp)
        os.system("phonopy --abinit -f %s" % (outs))

        # plot phonon band
        with open("band.conf", 'w') as fout:
            fout.write("ATOM_NAME =")
            for element in self.system.xyz.specie_labels:
                fout.write(" %s" % element)
            fout.write("\n")
            fout.write("DIM = %s\n" % supercell_n)
            fout.write("BAND = 0.5 0.5 0.5 0.0 0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.0\n")
        os.system("phonopy --abinit -c %s -p band.conf" % (pos_inp_name))
    #
    
    def get_disps(self, directory="./"):
        os.chdir(directory)
        os.system("ls | grep 'supercell-' > geo.data")
        disps = []
        with open("geo.data", 'r') as fin:
            for line in fin:
                disps.append(line.split(".")[0].split("-")[1])
        return disps
