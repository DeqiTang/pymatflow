#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg
from base.xyz import cp2k_xyz

"""
Usage:
    python converge_cutoff_cp2k.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.

    目前此脚本只进行pdos计算, 还没有实现结果分析
    cp2k官网是有Python脚本可以从 file.pdos 文件中提取数据的
    但是我看了下是用python2写的, 而且也很简单, 后面打算对其
    进行重写.
    http://wiki.wpi.edu/deskinsgroup/Density_of_States
    get-smearing-pdos.py
    pdos.py
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

base_project_name = "Calculate-PDOS"

if os.path.exists("./tmp-pdos"):
    shutil.rmtree("./tmp-pdos")
os.mkdir("./tmp-pdos")
os.chdir("./tmp-pdos")
shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")

inp_name = "pdos-calc.inp"

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
    fout.write("\t\t&PRINT\n")
    fout.write("\t\t\t&PDOS\n")
    fout.write("\t\t\t\tNLUMO -1\n") # print all projected DOS available
    fout.write("\t\t\t\tCOMPONENTS\n") # split the density by quantum number
    fout.write("\t\t\t&END PDOS\n")
    fout.write("\t\t&END PRINT\n")
    fout.write("\t&END DFT\n")
    # end dft
    fout.write("&END FORCE_EVAL\n")

# run the simulation
out_f_name = "pdos-calc.out"
os.system("cp2k.psmp -in %s > %s" % (inp_name, out_f_name))


# analyse the result

import matplotlib.pyplot as plt
