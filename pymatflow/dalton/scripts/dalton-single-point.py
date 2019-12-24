#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from emuhelper.dalton.base.xyz import dalton_xyz

"""
Usage:
    python single_point_dalton.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.
"""



# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]


xyz = dalton_xyz(sys.argv[1])

#base_project_name = "test"

if os.path.exists("./tmp-single-point"):
    shutil.rmtree("./tmp-single-point")
os.mkdir("./tmp-single-point")
os.chdir("./tmp-single-point")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])


comment_1 = "xxx"
comment_2 = "xxx"

mol_name = "single-point.mol"
dal_name = "single-point.dal"
with open(mol_name, 'w') as fout:
    fout.write("BASIS\n")
    fout.write("6-31G**\n")
    fout.write("%s\n" % comment_1)
    fout.write("%s\n" % comment_2)
xyz.to_dalton(mol_name)

with open(dal_name, 'w') as fout:
    fout.write("**DALTON INPUT\n")
    fout.write(".RUN PROPERTIES\n")
    fout.write("**WAVE FUNCTIONS\n")
    fout.write(".HF\n")
    fout.write(".MP2\n")
    fout.write("*SCF INPUT\n")
    fout.write(".DOUBLY OCCUPIED\n")
    fout.write(" 3 1 1 0\n")
    fout.write("*CONFIGURATION INPUT\n")
    fout.write("**START\n")
    fout.write("**PROPERTIES\n")
    fout.write(".DIPGRA\n")
    fout.write(".VIBANA\n")
    fout.write("**END OF DALTON INPUT\n")

# run the simulation
#os.system("dalton -mol single-point -dal single-point")
os.system("dalton -mol %s -dal %s" % (mol_name, dal_name))
