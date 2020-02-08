#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

from emuhelper.siesta.base.xyz import siesta_xyz


"""
Usage:
    python converge_ecut.py xxx.xyz ecut_min ecut_max ecut_step
    xxx.xyz is the input structure file

    make sure the xyz structure file and the corresponding pseudopotential
    file for the elements of the system is in the directory.
"""



# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
ecut_max = int(sys.argv[3])
ecut_step = int(sys.argv[4])


xyz = siesta_xyz(sys.argv[1])

system_name = "Test Meshcutoff Value"
system_label = "TestEcut"


if os.path.exists("./tmp-ecut"):
    shutil.rmtree("./tmp-ecut")
os.mkdir("./tmp-ecut")
os.chdir("./tmp-ecut")
#shutil.copyfile("../H.psf", "H.psf")
#shutil.copyfile("../Li.psf", "Li.psf")
#shutil.copyfile("../C.psf", "C.psf")
#shutil.copyfile("../Na.psf", "Na.psf")
#shutil.copyfile("../K.psf", "K.psf")
#shutil.copyfile("../Nb.psf", "Nb.psf")
#shutil.copyfile("../O.psf", "O.psf")
os.system("cp ../*.psf ./")


n_test = int((ecut_max - ecut_min) / ecut_step)
for i in range(n_test + 1):
    meshcutoff = int(ecut_min + i * ecut_step)
    fdf_name = "test-ecut-%d.fdf" % meshcutoff
    with open(fdf_name, 'w') as fout:
        fout.write("SystemName\t%s\n" % (system_name + str(meshcutoff)))
        fout.write("SystemLabel\t%s\n" % (system_label + str(meshcutoff)))
        fout.write("NumberOfSpecies\t%d\n" %  xyz.nspecies)
        fout.write("NumberOfAtoms\t%d\n" % xyz.natom)        
        fout.write("\n")
    xyz.to_fdf(fdf_name)
    #---------------------------------
    # DFT settings:
    # --------------------------------
    with open(fdf_name, 'a') as fout:
        fout.write("# DFT settings\n")
        fout.write("\n")
        fout.write("XC.functional  GGA\n")
        fout.write("XC.authors  PBE\n")
        fout.write("DM.Tolerance  1.d-6\n")
        fout.write("MaxSCFIterations  300\n")
        fout.write("DM.MixingWeight  0.1\n")
        fout.write("DM.NumberPulay  5\n")
        fout.write("DM.AllowExtrapolation  true\n")
        # do not DM.UseSaveDM so that different ecut does not affect each other
        fout.write("DM.UseSaveDM  false\n")         
        fout.write("SolutionMethod diagon\n")
        fout.write("MeshCutoff  %d Ry\n" % meshcutoff)
        fout.write("\n")

# run the simulation
for i in range(n_test + 1):
    meshcutoff = int(ecut_min + i * ecut_step)
    fdf_name = "test-ecut-%d.fdf" % meshcutoff
    out_f_name = "test-ecut-%d.out" % meshcutoff
    os.system("siesta < %s | tee %s" % (fdf_name, out_f_name))


# analyse the result
for i in range(n_test + 1):
    meshcutoff = int(ecut_min + i * ecut_step)
    out_f_name = "test-ecut-%d.out" % meshcutoff
    os.system("cat %s | grep 'Total =' >> energy-ecut.data" % out_f_name)

ecut = [ ecut_min + i * ecut_step for i in range(n_test + 1)]
energy = []
with open("energy-ecut.data", 'r') as fin:
    for line in fin:
        energy.append(float(line.split()[3]))

import matplotlib.pyplot as plt

plt.plot(ecut, energy)
plt.show()
