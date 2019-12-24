#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from emuhelper.orca.base.xyz import orca_xyz

"""
Usage:
    python converge_multiplicity_orca.py xxx.xyz n_test
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.
Note:
    设置自旋多重度有几个方式:
    1. 依据经验
    2. 依据晶体场理论
    3. 测试不同的自旋多重度值, 找出能量最低的那个
    对于普通的有机分子等无单电子体系, 是闭壳层, 自旋多重度为1
    可以验证此时计算得出能量比其它自旋多重度都要小, 其它自旋多重度对应激发态.
    对于含有过渡金属的体系, 自旋多重度的确定就比较复杂了, 这个脚本通过测试
    不同自旋多重度对应能量的值来找到合适的自旋多重度.
    规则是:
        电子数为奇数时测试自旋多重度为2, 4, 6...,
        电子数为偶数时测试自旋多重度为1, 3, 5...,
"""



# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]


xyz = orca_xyz(sys.argv[1])

#base_project_name = "test"

if os.path.exists("./tmp-converge-multiplicity"):
    shutil.rmtree("./tmp-converge-multiplicity")
os.mkdir("./tmp-converge-multiplicity")
os.chdir("./tmp-converge-multiplicity")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])


n_test = int(sys.argv[2])
n_electrons = 0
for atom in xyz.atoms:
    n_electrons += mg.Element(atom.name).number

if n_electrons % 2 == 1:
    multiplicity = [2*i+2 for i in range(n_test)]
else:
    multiplicity = [2*i+1 for i in range(n_test)]


for i in range(n_test):
    inp_name = "multiplicity-%d.inp" % multiplicity[i]
    with open(inp_name, 'w') as fout:
        fout.write("# comments\n")
        #fout.write("! def2-TZVP/C B3LYP\n")
        fout.write("! def2-SVP B3LYP\n")
        fout.write("! Angs\n")
        fout.write("* xyz 0 %d\n" % multiplicity[i])
    xyz.to_orca(inp_name)
    with open(inp_name, 'a') as fout:
        fout.write("*\n")

# run the simulation
for i in range(n_test):
    inp_name = "multiplicity-%d.inp" % multiplicity[i]
    os.system("orca %s | tee log-%d.txt" % (inp_name, multiplicity[i]))

# analysis result
import matplotlib.pyplot as plt
for i in range(n_test):
    os.system("cat log-%d.txt | grep 'Total Energy' >> multiplicity.data" % multiplicity[i])

energies = []
with open("multiplicity.data", 'r') as fin:
    for line in fin:
        energies.append(line.split()[5])

plt.plot(multiplicity, energies)
plt.show()
