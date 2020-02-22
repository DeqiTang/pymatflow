#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

from emuhelper.orca.base.xyz import orca_xyz

"""
Usage:
    python single_point_orca.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.
"""



# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]


xyz = orca_xyz(sys.argv[1])

#base_project_name = "test"

if os.path.exists("./tmp-single-point"):
    shutil.rmtree("./tmp-single-point")
os.mkdir("./tmp-single-point")
os.chdir("./tmp-single-point")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])


inp_name = "single-point.inp"
with open(inp_name, 'w') as fout:
    fout.write("# comments\n")
    fout.write("! def2-TZVP/C B3LYP\n")
    fout.write("! Angs\n")
    fout.write("* xyz 0 1\n")
xyz.to_orca(inp_name)
with open(inp_name, 'a') as fout:
    fout.write("*\n")
# run the simulation
os.system("orca %s | tee log.txt" % inp_name)
