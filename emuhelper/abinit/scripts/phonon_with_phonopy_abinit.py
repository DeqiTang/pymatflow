#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from emuhelper.abinit.base.xyz import abinit_xyz


"""
Usage:
    python phonon_with_phonopy_abinit.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential file
    for all the elements of the system is in the directory.
Note:
    参考: https://atztogo.github.io/phonopy/abinit.html
"""

        

# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

#ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#ecut_max = int(sys.argv[3])
#ecut_step = int(sys.argv[4])
cutoff = 40

supercell_n = "1 1 1"

xyz = abinit_xyz(sys.argv[1])

base_project_name = "phonon-calc"
if os.path.exists("./tmp-phonon-with-phonopy"):
    shutil.rmtree("./tmp-phonon-with-phonopy")
os.mkdir("./tmp-phonon-with-phonopy")
os.chdir("./tmp-phonon-with-phonopy")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")
#shutil.copyfile("../Li.psp8", "Li.psp8")
#shutil.copyfile("../H.psp8", "H.psp8")
os.system("cp ../*.psp8 ./")

head_inp_name = "head-phonon.in"
with open(head_inp_name, 'w') as fout:
    fout.write("kptopt 1\n")
    fout.write("ngkpt 1 1 1\n")
    fout.write("occopt 3\n")
    fout.write("ecut %d\n" % cutoff)
    fout.write("toldfe 1.0d-6\n")
    fout.write("nstep 100\n")
    fout.write("diemac 2.0\n")
    fout.write("\n")

pos_inp_name = "pos.in"
xyz.to_abinit(pos_inp_name)
os.system("phonopy --abinit -d --dim='%s' -c %s" % (supercell_n, pos_inp_name))

os.system("ls | grep 'supercell-' > geo.data")
disp_dirs = []
with open("geo.data", 'r') as fin:
    for line in fin:
        disp_dirs.append(line.split(".")[0].split("-")[1])

for disp in disp_dirs:
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
    os.system("head %s -n -%d >> ./disp-%s/supercell-%s.in" % (pos_inp_name, (xyz.natom + 2), disp, disp))
    os.system("tail supercell-%s.in -n %d >> ./disp-%s/supercell-%s.in" % (disp, (xyz.natom + 1), disp, disp))
    with open("./disp-%s/supercell-%s.files" % (disp, disp), 'w') as fout:
        fout.write("supercell-%s.in\n" % disp)
        fout.write("supercell-%s.out\n" % disp)
        fout.write("supercell-%si\n" % disp)
        fout.write("supercell-%so\n" % disp)
        fout.write("temp\n")
        for element in xyz.specie_labels:
            fout.write("%s\n" % (element + ".psp8"))

# run the simulation
for disp in disp_dirs:
    os.chdir("disp-%s" % disp)
    os.system("abinit < supercell-%s.files" % disp)
    os.chdir("../")

# analyse the result

import matplotlib.pyplot as plt

os.system("phonopy --abinit -f disp-{001..%s}/supercell-*.out" % (disp_dirs[-1]))

# plot phonon band
with open("band.conf", 'w') as fout:
    fout.write("ATOM_NAME =")
    for element in xyz.specie_labels:
        fout.write(" %s" % element)
    fout.write("\n")
    fout.write("DIM = %s\n" % supercell_n)
    fout.write("BAND = 0.5 0.5 0.5 0.0 0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.0\n")
os.system("phonopy --abinit -c %s -p band.conf" % (pos_inp_name))
