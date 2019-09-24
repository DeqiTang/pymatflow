#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from base.xyz import abinit_xyz

"""
Usage:
    python converge_ecut.py xxx.xyz ecut_min ecut_max ecut_step
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential file
    for all the elements of the system is in the directory.
"""



# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
ecut_max = int(sys.argv[3])
ecut_step = int(sys.argv[4])


xyz = abinit_xyz(sys.argv[1])

base_project_name = "test"
if os.path.exists("./tmp-ecut"):
    shutil.rmtree("./tmp-ecut")
os.mkdir("./tmp-ecut")
os.chdir("./tmp-ecut")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")
#shutil.copyfile("../Li.psp8", "Li.psp8")
#shutil.copyfile("../H.psp8", "H.psp8")
os.system("cp ../*.psp8 ./")


n_test = int((ecut_max - ecut_min) / ecut_step)
for i in range(n_test + 1):
    cutoff = int(ecut_min + i * ecut_step)
    rel_cutoff = cutoff / 3
    inp_name = "test-ecut-%d.in" % cutoff
    files_name = "test-ecut-%d.files" % cutoff
    with open(files_name, 'w') as fout:
        fout.write(inp_name)
        fout.write("\n")
        fout.write("test-ecut-%d.out\n" % cutoff)
        fout.write("test-ecut-%di\n" % cutoff)
        fout.write("test-ecut-%do\n" % cutoff)
        fout.write("temp\n")
        for element in xyz.specie_labels:
            fout.write("%s\n" % (element + ".psp8"))
        #
    with open(inp_name, 'w') as fout:
        fout.write("ecut %d\n" % cutoff)
        fout.write("ixc = 11\n")
        fout.write("kptopt 1\n")
        fout.write("ngkpt 1 1 1\n")
        fout.write("nstep 100\n")
        fout.write("toldfe 1.0d-6\n")
        fout.write("diemac 2.0\n")
        fout.write("\n")
    xyz.to_abinit(inp_name)
# run the simulation
for i in range(n_test + 1):
    cutoff = int(ecut_min + i * ecut_step)
    files_name = "test-ecut-%d.files" % cutoff
    os.system("abinit < %s" % (files_name))

# analyse the result
for i in range(n_test + 1):
    cutoff = int(ecut_min + i * ecut_step)
    out_f_name = "test-ecut-%d.out" % cutoff
    os.system("cat %s | grep 'Etotal=' >> energy-ecut.data" % out_f_name)

ecut = [ ecut_min + i * ecut_step for i in range(n_test + 1)]
energy = []
with open("energy-ecut.data", 'r') as fin:
    for line in fin:
        energy.append(float(line.split()[2]))

import matplotlib.pyplot as plt

#for i in range(len(energy)):
#    energy[i] = energy[i] - 31
plt.plot(ecut, energy)
plt.show()
