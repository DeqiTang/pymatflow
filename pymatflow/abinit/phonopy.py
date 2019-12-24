#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import seekpath
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.abinit.base.electrons import abinit_electrons
from pymatflow.abinit.base.system import abinit_system

"""
Usage:
Note:
    参考: https://atztogo.github.io/phonopy/abinit.html
"""

        
class phonopy_run:
    """
    """
    def __init__(self, xyz_f):
        self.system = abinit_system(xyz_f)
        self.electrons = abinit_electrons()
    
        self.electrons.basic_setting()

        self.supercell_n = [1, 1, 1]

    def phonopy(self, directory="tmp-abinit-phonopy", head_inpname="head-phonon.in", pos_inpname="pos.in", mpi="", runopt="gen",
            electrons={}, kpoints={}, supercell_n=[1, 1, 1]):
        """
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
            os.system("cp %s %s/" % (self.system.xyz.file, directory))
           
            self.supercell_n = supercell_n

            self.electrons.set_scf_nscf("scf")
            self.electrons.set_params(electrons)
            self.electrons.kpoints.set_params(kpoints)
            #
            with open(os.path.join(directory, head_inpname), 'w') as fout:
                self.electrons.to_in(fout)
            with open(os.path.join(directory, pos_inpname), 'w') as fout:
                self.system.to_in(fout)
   
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
                os.system("head %s -n -%d >> ./disp-%s/supercell-%s.in" % (pos_inpname, (self.system.xyz.natom + 2), disp, disp))
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

            # generate yhbatch scripts
            with open(os.path.join(directory, "phonopy-job.sub"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                for disp in disps:
                    fout.write("cd disp-%s\n" % disp)
                    fout.write("yhrun -N 1 -n 24 abinit < supercell-%s.files\n" % disp)
                    fout.write("cd ../\n")
            
            # generate the result analysis bash scripts and necessary config files
            os.chdir(directory)
                        
            with open("mesh.conf", 'w') as fout:
                fout.write("ATOM_NAME =")
                for element in self.system.xyz.specie_labels:
                    fout.write(" %s" % element)
                fout.write("\n")
                fout.write("DIM = %d %d %d\n" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))
                fout.write("MP = 8 8 8\n")
            
                    
            with open("pdos.conf", 'w') as fout:
                fout.write("ATOM_NAME =")
                for element in self.system.xyz.specie_labels:
                    fout.write(" %s" % element)
                fout.write("\n")
                fout.write("DIM = %d %d %d\n" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))
                fout.write("MP = 8 8 8\n")
                fout.write("PDOS = 1 2, 3 4 5 5\n")
            
            
            with open("band.conf", 'w') as fout:
                fout.write("ATOM_NAME =")
                for element in self.system.xyz.specie_labels:
                    fout.write(" %s" % element)
                fout.write("\n")
                # the use of PRIMITIVE_AXES will find the primitive cell of the structure
                # and use it to analyse the phonon band structure
                # however, the use of primitive cell will not affect the q path setting
                # so whether we use PRIMITIVE cell or not, we can set the same q path
                fout.write("PRIMITIVE_AXES = AUTO\n") # we can also specify a matrix, but AUTO is recommended now in phonopy
                fout.write("GAMMA_CENTER = .TRUE.\n")
                fout.write("BAND_POINTS = 101\n")
                fout.write("BAND_CONNECTION = .TRUE.\n")

                fout.write("DIM = %d %d %d\n" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))
                #fout.write("BAND = 0.5 0.5 0.5 0.0 0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.0\n")
                fout.write("BAND =")
                # --------------
                # using seekpath to set q path
                # --------------
                lattice = self.system.xyz.cell  # [self.system.xyz.cell[0:3], self.system.xyz.cell[3:6], self.system.xyz.cell[6:9]]
                positions = []
                numbers = []
                #a = np.sqrt(self.system.xyz.cell[0]**2 + self.system.xyz.cell[1]**2 + self.system.xyz.cell[2]**2)
                #b = np.sqrt(self.system.xyz.cell[3]**2 + self.system.xyz.cell[4]**2 + self.system.xyz.cell[5]**2)
                #c = np.sqrt(self.system.xyz.cell[6]**2 + self.system.xyz.cell[7]**2 + self.system.xyz.cell[8]**2)
                a = np.sqrt(self.system.xyz.cell[0][0]**2 + self.system.xyz.cell[0][1]**2 + self.system.xyz.cell[0][2]**2)
                b = np.sqrt(self.system.xyz.cell[1][0]**2 + self.system.xyz.cell[1][1]**2 + self.system.xyz.cell[1][2]**2)
                c = np.sqrt(self.system.xyz.cell[2][0]**2 + self.system.xyz.cell[2][1]**2 + self.system.xyz.cell[2][2]**2)
                for atom in self.system.xyz.atoms:
                    positions.append([atom.x / a, atom.y / b, atom.z / c])
                    numbers.append(self.system.xyz.specie_labels[atom.name])
                structure = (lattice, positions, numbers)
                kpoints_seekpath = seekpath.get_path(structure)
                point = kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][0]]
                fout.write(" %f %f %f" % (point[0], point[1], point[2])) #self.arts..kpoints_seekpath["path"][0][0]
                point = kpoints_seekpath["point_coords"][kpoints_seekpath["path"][0][1]]
                fout.write(" %f %f %f" % (point[0], point[1], point[2])) #self.arts.kpoints_seekpath["path"][0][1]
                for i in range(1, len(kpoints_seekpath["path"])):
                    if kpoints_seekpath["path"][i][0] == kpoints_seekpath["path"][i-1][1]:
                        point = kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]]
                        fout.write(" %f %f %f" % (point[0], point[1], point[2])) #self.arts.kpoints_seekpath["path"][i][1]))
                    else:
                        point = kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][0]]
                        fout.write(" %f %f %f" % (point[0], point[1], point[2])) #self.arts.kpoints_seekpath["path"][i][0]))
                        point = kpoints_seekpath["point_coords"][kpoints_seekpath["path"][i][1]]
                        fout.write(" %f %f %f" % (point[0], point[1], point[2]))  #self.kpoints_seekpath["path"][i][1]))
                fout.write("\n")
                fout.write("BAND_LABELS =")
                point = kpoints_seekpath["path"][0][0]
                if point == "GAMMA":
                    fout.write(" $\Gamma$")
                else:
                    fout.write(" $%s$" % point)
                point = kpoints_seekpath["path"][0][1]
                if point == "GAMMA":
                    fout.write(" $\Gamma$")
                else:
                    fout.write(" $%s$" % point)
                for i in range(1, len(kpoints_seekpath["path"])):
                    if kpoints_seekpath["path"][i][0] == kpoints_seekpath["path"][i-1][1]:
                        point = kpoints_seekpath["path"][i][1]
                        if point == "GAMMA":
                            fout.write(" $\Gamma$")
                        else:
                            fout.write(" $%s$" % point)
                    else:
                        point = kpoints_seekpath["path"][i][0]
                        if point == "GAMMA":
                            fout.write(" $\Gamma$")
                        else:
                            fout.write(" $%s$" % point)
                        point = kpoints_seekpath["path"][i][1]
                        if point == "GAMMA":
                            fout.write(" $\Gamma$")
                        else:
                            fout.write(" $%s$" % point)
                fout.write("\n")
                #
            #
            with open("phonopy-analysis.sh", 'w') as fout:
                fout.write("#!/bin/bash\n\n")
                outs = ""
                for disp in disps:
                    outs += "disp-%s/supercell-%s.out " % (disp, disp)
                fout.write("# generate FORCE_SETS\n")
                fout.write("phonopy --abinit -f %s\n" % (outs))
                fout.write("# plot the density of states (DOS)\n")
                fout.write("phonopy --abinit -p mesh.conf -c %s\n" % pos_inpname)
                fout.write("# Thermal properties are calculated with the sampling mesh by:\n")
                fout.write("phonopy --abinit -t mesh.conf -c %s\n" % pos_inpname)
                fout.write("# Thermal properties can be plotted by:\n")
                fout.write("phonopy --abinit -t -p mesh.conf -c %s\n" % pos_inpname)
                fout.write("# calculate Projected DOS and plot it\n")
                fout.write("phonopy --abinit -p pdos.conf -c %s\n" % pos_inpname)
                fout.write("# plot phonon band\n")
                fout.write("phonopy --abinit -c %s -p band.conf\n" % (pos_inpname))
            os.chdir("../")
            # end generate the result analysis bash script and necessary config giles

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
