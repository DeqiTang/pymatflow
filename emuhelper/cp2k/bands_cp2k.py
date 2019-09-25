#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from emuhelper.cp2k.base.xyz import cp2k_xyz

"""
Usage:
    python bands_cp2k.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.

    能带的计算需要设置K点采样
    参考SeeK-path这个工具:https://tools.materialscloud.org/seekpath/
    可以输入结构自动为cp2k, qe等程序生成k点

    此脚本会使cp2k生成xxx.bs能带文件
    可以用cp2k官网出现的一个cp2k_bs2csv.py进行处理后用gnuplot绘图
    我考虑后期自行实现该处理, 并使用Python作图
    参考: https://www.cp2k.org/exercises:2018_uzh_cmest:pdos
"""


        

# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

#cutoff_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#cutoff_max = int(sys.argv[3])
#cutoff_step = int(sys.argv[4])
#rel_cutoff = int(sys.argv[5])
cutoff = 60 
rel_cutoff = 30

xyz = cp2k_xyz(sys.argv[1])

base_project_name = "Calculate-Bands"

if os.path.exists("./tmp-bands"):
    shutil.rmtree("./tmp-bands")
os.mkdir("./tmp-bands")
os.chdir("./tmp-bands")
shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")

inp_name = "bands-calc.inp"

with open(inp_name, 'w') as fout:
    fout.write("&GLOBAL\n")
    fout.write("\tPROJECT\t%s\n" % (base_project_name))
    fout.write("\tRUN_TYPE ENERGY\n")
    fout.write("\tPRINT_LEVEL MEDIUM\n")
    fout.write("&END GLOBAL\n")
    fout.write("\n")

    fout.write("&FORCE_EVAL\n")
    fout.write("\tMETHOD Quickstep\n")
# subsys
xyz.to_subsys(inp_name)
# end subsys
with open(inp_name, 'a') as fout:
    # dft
    fout.write("\t&DFT\n")
    fout.write("\t\tBASIS_SET_FILE_NAME BASIS_SET\n")
    fout.write("\t\tPOTENTIAL_FILE_NAME GTH_POTENTIALS\n")
    fout.write("\t\t&QS\n")
    fout.write("\t\t\tEPS_DEFAULT 1.0E-10\n")
    fout.write("\t\t\tEXTRAPOLATION USE_GUESS\n") # required for K-Point sampling
    fout.write("\t\t&END QS\n")
    fout.write("\t\t&MGRID\n")
    fout.write("\t\t\tNGRIDS 4\n")
    fout.write("\t\t\tCUTOFF %d\n" % cutoff)
    fout.write("\t\t\tREL_CUTOFF %d\n" % rel_cutoff)
    fout.write("\t\t&END MGRID\n")
    fout.write("\t\t&XC\n")
    fout.write("\t\t\t&XC_FUNCTIONAL PADE\n")
    fout.write("\t\t\t&END XC_FUNCTIONAL\n")
    fout.write("\t\t&END XC\n")
    fout.write("\t\t&SCF\n")
    fout.write("\t\t\tSCF_GUESS ATOMIC\n")
    fout.write("\t\t\tEPS_SCF 1.0E-06\n")
    fout.write("\t\t\tMAX_SCF 300\n")
    fout.write("\t\t\t&DIAGONALIZATION ON\n")
    fout.write("\t\t\t\tALGORITHM STANDARD\n")
    fout.write("\t\t\t&END DIAGONALIZATION\n")
    fout.write("\t\t\t&MIXING T\n")
    fout.write("\t\t\t\tMETHOD BROYDEN_MIXING\n")
    fout.write("\t\t\t\tALPHA 0.4\n")
    fout.write("\t\t\t\tNBROYDEN 8\n")
    fout.write("\t\t\t&END MIXING\n")
    fout.write("\t\t&END SCF\n")

    fout.write("\t\t&POISSON\n")
    fout.write("\t\t\tPERIODIC XYZ\n")
    fout.write("\t\t&END POISSON\n")

    fout.write("\t\t&PRINT\n")
    fout.write("\t\t\t&BAND_STRUCTURE\n")
    fout.write("\t\t\t\tADDED_MOS 2\n") 
    fout.write("\t\t\t\tFILE_NAME bands.bs\n")
    fout.write("\t\t\t\t&KPOINT_SET\n")
    fout.write("\t\t\t\t\tUNITS B_VECTOR\n")
    fout.write("\t\t\t\t\tSPECIAL_POINT 0.0000000000     0.0000000000     0.0000000000 #GAMA\n")
    fout.write("\t\t\t\t\tSPECIAL_POINT 0.5000000000     0.0000000000     0.5000000000 #X\n")
    fout.write("\t\t\t\t\tSPECIAL_POINT 0.5000000000     0.2500000000     0.7500000000 #W\n")
    fout.write("\t\t\t\t\tSPECIAL_POINT 0.0000000000     0.0000000000     0.0000000000 #GAMA\n")
    fout.write("\t\t\t\t\tSPECIAL_POINT 0.6250000000     0.2500000000     0.6250000000 #U\n")
    fout.write("\t\t\t\t\tSPECIAL_POINT 0.5000000000     0.2500000000     0.7500000000 #W\n")
    fout.write("\t\t\t\t\tNPOINTS 10\n")
    fout.write("\t\t\t\t&END KPOINT_SET\n")
    fout.write("\t\t\t&END BAND_STRUCTURE\n")
    fout.write("\t\t&END PRINT\n")
    fout.write("\t&END DFT\n")
    # end dft
    fout.write("&END FORCE_EVAL\n")

# run the simulation
out_f_name = "bands-calc.out"
os.system("cp2k.psmp -in %s > %s" % (inp_name, out_f_name))


# analyse the result

import matplotlib.pyplot as plt
