#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

from emuhelper.abinit.base.xyz import abinit_xyz

"""
Usage:
    python bands_abinit.py xxx.xyz
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

base_project_name = "bands-calc"
if os.path.exists("./tmp-bands"):
    shutil.rmtree("./tmp-bands")
os.mkdir("./tmp-bands")
os.chdir("./tmp-bands")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")
#shutil.copyfile("../Li.psp8", "Li.psp8")
#shutil.copyfile("../H.psp8", "H.psp8")
os.system("cp ../*.psp8 ./")

inp_name = "bands-calc.in"
files_name = "bands-calc.files"
with open(files_name, 'w') as fout:
    fout.write(inp_name)
    fout.write("\n")
    fout.write("bands-calc.out\n")
    fout.write("bands-calci\n")
    fout.write("bands-calco\n")
    fout.write("temp\n")
    for element in xyz.specie_labels:
        fout.write("%s\n" % (element + ".psp8"))
    #
with open(inp_name, 'w') as fout:
    fout.write("ndtset 2\n")
    # dataset 1: scf calculation
    fout.write("# Dataset 1: scf\n")
    fout.write("kptopt1 1\n")
    fout.write("nshiftk1 4\n")
    fout.write("shiftk1 0.5 0.5 0.5\n")
    fout.write("  0.5 0.0 0.0\n")
    fout.write("  0.0 0.5 0.0\n")
    fout.write("  0.0 0.0 0.5\n")
    fout.write("ngkpt1 4 4 4\n")
    fout.write("prtden1 1\n")
    fout.write("toldfe1 1.0d-6\n")
    # dataset 2: band calculation
    fout.write("iscf2 -2\n")
    fout.write("getden2 -1\n")
    fout.write("kptopt2 -3\n")
    fout.write("nband2 8\n")
    fout.write("ndivk2 10 12 17\n")
    fout.write("kptbounds2 0.5 0.0 0.0 # L point\n")
    fout.write("0.0 0.0 0.0 # Gamma\n")
    fout.write("0.0 0.5 0.5 # X\n")
    fout.write("1.0 1.0 1.0 # Gamma in another cell\n")
    fout.write("tolwfr2 1.0d-12\n")
    fout.write("enunit2 1\n")
    fout.write("\n")

    fout.write("ecut %d\n" % cutoff)
    fout.write("occopt 3\n")
    #fout.write("prtdos 1\n") # prtdos = 1, 2, 3
    fout.write("nstep 100\n")
    fout.write("diemac 2.0\n")
    fout.write("\n")
xyz.to_abinit(inp_name)

# run the simulation
os.system("abinit < %s" % (files_name))

# analyse the result

import matplotlib.pyplot as plt
os.system("abiopen.py bands-calco_DS1_GSR.nc -e -sns=talk")
