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
    python geo_opt_orca.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.
"""


        

# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]


xyz = orca_xyz(sys.argv[1]) 

#base_project_name = "test"

if os.path.exists("./tmp-geo-opt"):
    shutil.rmtree("./tmp-geo-opt")
os.mkdir("./tmp-geo-opt")
os.chdir("./tmp-geo-opt")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])


inp_name = "geo-opt.inp"
with open(inp_name, 'w') as fout:
    fout.write("# comments\n")
    fout.write("! def2-TZVP B3LYP/G Opt TightOpt\n")
    # 注意自旋多重度的正确设置很重要:对计算时间和结果的影响都非常巨大
    fout.write("! Angs\n")
    fout.write("* xyz 0 1\n") 
xyz.to_orca(inp_name)
with open(inp_name, 'a') as fout:
    fout.write("*\n")
# run the simulation
out_name = "geo-opt.out"
os.system("orca %s | tee %s" % (inp_name, out_name))
