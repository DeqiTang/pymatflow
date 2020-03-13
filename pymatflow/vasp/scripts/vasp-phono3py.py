#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil


from pymatflow.vasp.base.xyz import vasp_studio

"""
Usage:
    python phonon3py_vasp.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and the POTCAR is in the directory.
    make sure the element in the POTCAR is in order of increasing the atom number of
    the element: 确保POTCAR中元素排列的顺序是原子序数从小到大
Note:
    参考: https://atztogo.github.io/phono3py/
    pip install --user phono3py # 可能需要 sudo pacman -S --noconfirm lapacke
    目前脚本已经编写好, 但是无法测试, 因为phono3py似乎还有问题第一步生成
    生成有位移的POSCAR的时候就会出错. 过段时间再进行测试
"""


# =======================================
# Constructing the input file for VASP
# =======================================

cutoff = 125

supercell_n = "1 1 1"

xyz = vasp_xyz(sys.argv[1])

base_project_name = "Phonon Calculation"

if os.path.exists("./tmp-phono3py"):
    shutil.rmtree("./tmp-phono3py")
os.mkdir("./tmp-phono3py")
os.chdir("./tmp-phono3py")

# now I generate POTCAR manually, but actually I can use the command
# tool pmg from pymatgen to automatically generate POTCAR
shutil.copyfile("../POTCAR", "POTCAR")

xyz.to_poscar("POSCAR")

# ======================================
# Phonopy set up
# ======================================

# constructing the INCAR for the phonon calculation
with open("INCAR", 'w') as fout:
    fout.write("SYSTEM = XXX\n")
    fout.write("ISMEAR = 0\n") # 此设置可能经常需要修改
    fout.write("SIGMA = 0.2\n")
    fout.write("IBRION = -1\n")
    fout.write("LREAL = .FALSE.\n")
    fout.write("LWAVE = .FALSE.\n")
    fout.write("LCHARG = .FALSE.\n")
    fout.write("ENCUT = %d\n" % cutoff)
    fout.write("\n")
#
with open("KPOINTS", 'w') as fout:
    fout.write("K-POINTS\n")
    fout.write("0\n") # automatically generate the grid
    fout.write("Gamma\n") # Gamma centered
    fout.write("1 1 1\n")
    fout.write("0 0 0\n") # usually set them to 0

## Construct and run every POSCAR scf
os.system("phono3py -d --dim='%s' -c POSCAR" % supercell_n)
os.system("ls | grep 'POSCAR-' > pos.data")
disp_dirs = []
with open("pos.data", 'r') as fin:
    for line in fin:
        disp_dirs.append(line.split("\n")[0].split("-")[1])

for disp in disp_dirs:
    os.mkdir("disp-%s" % (disp))
    os.chdir("disp-%s" % (disp))
    shutil.copyfile("../POSCAR-%s" % disp, "POSCAR")
    shutil.copyfile("../INCAR", "INCAR")
    shutil.copyfile("../POTCAR", "POTCAR")
    shutil.copyfile("../KPOINTS", "KPOINTS")
    os.chdir("../")

for disp in disp_dirs:
    os.chdir("disp-%s" % disp)
    os.system("vasp_std")
    os.chdir("../")


# analyse the result
os.system("phono3py --cf3 disp-{001..%s}/vasprun.xml" % (disp_dirs[-1]))
os.system("phono3py --dim='%s' -c POSCAR" % supercell_n)

# Thermal conductivity calculation
os.system("phono3py --fc3 --dim='%s' --mesh='11 11 11' -c POSCAR --br" % supercell_n)
