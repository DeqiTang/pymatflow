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
    python geo_opt_abinit.py xxx.xyz
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

xyz = abinit_xyz(sys.argv[1])

base_project_name = "geo-opt-calc"
if os.path.exists("./tmp-geo-opt"):
    shutil.rmtree("./tmp-geo-opt")
os.mkdir("./tmp-geo-opt")
os.chdir("./tmp-geo-opt")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")
#shutil.copyfile("../Li.psp8", "Li.psp8")
#shutil.copyfile("../H.psp8", "H.psp8")
os.system("cp ../*.psp8 ./")


inp_name = "geo-opt-calc.in"
files_name = "geo-opt-calc.files"
with open(files_name, 'w') as fout:
    fout.write(inp_name)
    fout.write("\n")
    fout.write("geo-opt-calc.out\n")
    fout.write("geo-opt-calci\n")
    fout.write("geo-opt-calco\n")
    fout.write("temp\n")
    for element in xyz.specie_labels:
        fout.write("%s\n" % (element + ".psp8"))
    #
with open(inp_name, 'w') as fout:
    fout.write("ecut %d\n" % cutoff)
    fout.write("ixc = 11\n")
    fout.write("kptopt 1\n")
    fout.write("ngkpt 1 1 1\n")
    fout.write("occopt 3\n") # fermi dirac smearing of occupation
    fout.write("nstep 100\n")
    fout.write("toldfe 1.0d-6\n")
    fout.write("diemac 2.0\n")
    # geo opt setting
    fout.write("ionmov 3\n")
    fout.write("optcell 0\n")
    fout.write("ntime 100\n")
    fout.write("tolmxde 0.0001 eV\n") # 
    fout.write("tolmxf 0 #5.0d-3  # Ha/Bohr\n") # tolmxde 与tolmxf两者只能设置一个为有效值, 其一需要定义为0
    fout.write("\n")
xyz.to_abinit(inp_name)

# run the simulation
#out_f_name = "geo-opt-calc.out.log"
os.system("abinit < %s" % (files_name))

# analyse the result

import matplotlib.pyplot as plt
# use abiview.py to view information in relaxation
os.system("abiview.py hist %so_HIST.nc" % base_project_name)


