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
    python dos_abinit.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential file
    for all the elements of the system is in the directory.
"""



# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

#ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#ecut_max = int(sys.argv[3])
#ecut_step = int(sys.argv[4])
cutoff = 40

scale_cart = "3*1" # '3*1' actually means no scaling

xyz = abinit_xyz(sys.argv[1])

base_project_name = "dos-calc"
if os.path.exists("./tmp-dos"):
    shutil.rmtree("./tmp-dos")
os.mkdir("./tmp-dos")
os.chdir("./tmp-dos")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")
#shutil.copyfile("../Li.psp8", "Li.psp8")
#shutil.copyfile("../H.psp8", "H.psp8")
os.system("cp ../*.psp8 ./")


inp_name = "dos-calc.in"
files_name = "dos-calc.files"
with open(files_name, 'w') as fout:
    fout.write(inp_name)
    fout.write("\n")
    fout.write("dos-calc.out\n")
    fout.write("dos-calci\n")
    fout.write("dos-calco\n")
    fout.write("temp\n")
    for element in xyz.specie_labels:
        fout.write("%s\n" % (element + ".psp8"))
    #
with open(inp_name, 'w') as fout:
    fout.write("scalecart %s\n" % scale_cart) # 一定要注意scalecart的重要性, 设置大了可能会耗尽内存, 但对于周期体系又和重要
    fout.write("ecut %d\n" % cutoff)
    # k-point setting
    fout.write("kptopt 1\n")
    fout.write("ngkpt 3 3 3\n")
    # dos setting
    fout.write("occopt 3\n")
    fout.write("prtdos 1\n") # prtdos = 1, 2, 3
    fout.write("nstep 100\n")
    fout.write("toldfe 1.0d-6\n")
    fout.write("diemac 2.0\n")
    fout.write("\n")
xyz.to_abinit(inp_name)

# run the simulation
os.system("abinit < %s" % (files_name))

# analyse the result
# import matplotlib.pyplot as plt
os.system("abiopen.py %s --expose -sns=talk" % (base_project_name + "o_GSR.nc"))

