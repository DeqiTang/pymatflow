#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import seekpath
import matplotlib.pyplot as plt

from pymatflow.abinit.abinit import abinit

"""
Usage:
Note:
    参考: https://atztogo.github.io/phonopy/abinit.html
"""


class phonopy_run(abinit):
    """
    """
    def __init__(self):
        super().__init__()

        #self.input.guard.set_queen(queen="static", electrons=self.input.electrons, system=self.input.system)

        self.input.electrons.basic_setting()

        self.supercell_n = [1, 1, 1]

    def phonopy(self, directory="tmp-abinit-phonopy", head_inpname="head-phonon.in", pos_inpname="pos.in", mpi="", runopt="gen",
        jobname="phonopy", nodes=1, ppn=32):
        """
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
            os.system("cp %s %s/" % (self.input.system.xyz.file, directory))


            self.input.electrons.set_scf_nscf("scf")
            #
            with open(os.path.join(directory, head_inpname), 'w') as fout:
                self.input.electrons.to_in(fout)
            with open(os.path.join(directory, pos_inpname), 'w') as fout:
                self.input.system.to_in(fout)

            os.chdir(directory)
            os.system("phonopy --abinit -d --dim='%d %d %d' -c %s" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2], pos_inpname))

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
                os.system("cat %s > ./disp-%s/supercell-%s.in" % (head_inpname, disp, disp))
                os.system("head %s -n -%d >> ./disp-%s/supercell-%s.in" % (pos_inpname, (self.input.system.xyz.natom + 2), disp, disp))
                os.system("tail supercell-%s.in -n %d >> ./disp-%s/supercell-%s.in" % (disp, (self.input.system.xyz.natom + 1), disp, disp))
                with open("./disp-%s/supercell-%s.files" % (disp, disp), 'w') as fout:
                    fout.write("supercell-%s.in\n" % disp)
                    fout.write("supercell-%s.out\n" % disp)
                    fout.write("supercell-%si\n" % disp)
                    fout.write("supercell-%so\n" % disp)
                    fout.write("temp\n")
                    for element in self.input.system.xyz.specie_labels:
                        fout.write("%s\n" % (element + ".psp8"))
            os.chdir("../")

            # generate yhbatch scripts
            with open(os.path.join(directory, "phonopy-job.sub"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                for disp in disps:
                    fout.write("cd disp-%s\n" % disp)
                    fout.write("yhrun -N 1 -n 24 abinit < supercell-%s.files\n" % disp)
                    fout.write("cd ../\n")

            # generate pbs scripts
            with open(os.path.join(directory, "phonopy-job.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % jobname)
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for disp in disps:
                    fout.write("cd disp-%s\n" % disp)
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE abinit < supercell-%s.files\n" % disp)
                    fout.write("cd ../\n")


        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            disps = self.get_disps("./")
            for disp in disps:
                os.chdir("disp-%s" % disp)
                os.system("abinit < supercell-%s.files" % disp)
                os.chdir("../")
            os.chdir("../")


    def get_disps(self, directory="./"):
        os.chdir(directory)
        os.system("ls | grep 'supercell-' > geo.data")
        disps = []
        with open("geo.data", 'r') as fin:
            for line in fin:
                disps.append(line.split(".")[0].split("-")[1])
        return disps
