#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import re
import pymatgen as mg

from emuhelper.qe.base.xyz import qe_xyz

"""
Usage:
    python converge_ecutwfc_pwscf.py xxx.xyz ecut_wfc_min ecut_wfc_max ecut_wfc_step
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential
    file for all the elements in the system is in the directory.
Note:
    here while finding the best ecutwfc I always set the the ecutrho to the default
    value: 4 * ecutwfc
"""



# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

ecut_wfc_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
ecut_wfc_max = int(sys.argv[3])
ecut_wfc_step = int(sys.argv[4])


xyz = qe_xyz(sys.argv[1])

title = "Test ecutwfc"
base_prefix = "knn"

if os.path.exists("./tmp-ecutwfc"):
    shutil.rmtree("./tmp-ecutwfc")
os.mkdir("./tmp-ecutwfc")
os.chdir("./tmp-ecutwfc")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")
os.system("cp ../*.UPF ./")

n_test = int((ecut_wfc_max - ecut_wfc_min) / ecut_wfc_step)
for i in range(n_test + 1):
    ecut_wfc = int(ecut_wfc_min + i * ecut_wfc_step)
    ecut_rho = ecut_wfc * 4 # using default value for ecut_rho: 4 * ecutwfc
    inp_name = "test-ecutwfc-%d.in" % ecut_wfc
    with open(inp_name, 'w') as fout:
        fout.write("&control\n")
        fout.write("calculation = 'scf'\n")
        fout.write("title = '%s'\n" % title)
        fout.write("prefix = '%s'\n" % (base_prefix + str(ecut_wfc)))
        fout.write("restart_mode = 'from_scratch'\n")
        fout.write("nstep = 300\n")
        fout.write("outdir = '%s'\n" % ("./tmp-" + str(ecut_wfc)))
        fout.write("pseudo_dir = './'\n")
        fout.write("wf_collect = .true.\n")
        fout.write("tstress = .true.\n")
        fout.write("tprnfor = .true.\n")
        fout.write("/\n")
        fout.write("\n")

        fout.write("&system\n")
        fout.write("ibrav = 0\n")
        fout.write("nat = %d\n" % xyz.natom)
        fout.write("ntyp = %d\n" % xyz.nspecies)
        fout.write("nspin = 1\n")
        #fout.write("nspin = 2\n")
        #fout.write("starting_magnetization(1) = 1\n")
        #fout.write("starting_magnetization(2) = 1\n")
        fout.write("ecutwfc = %d\n" % ecut_wfc)
        fout.write("ecutrho = %d \n" % ecut_rho)  # default value: 4 * ecutwfc
        fout.write("input_DFT = 'PBE'\n")
        fout.write("occupations = 'smearing'\n") # smearing, tetrahedra, fixed
        fout.write("degauss = 1.0d-4\n") # default: 0
        fout.write("smearing = 'gaussian'\n") # default is gaussian, and you can use fermi-dirac
        fout.write("/\n")
        fout.write("\n")

        fout.write("&electrons\n")
        fout.write("electron_maxstep = 300\n")
        #fout.write("conv_thr = 1.0d-10\n")
        fout.write("conv_thr = 1.0d-5\n")
        fout.write("mixing_mode = 'plain'\n")
        fout.write("mixing_beta = 0.3d0\n")
        fout.write("scf_must_converge = .true.\n")
        fout.write("/\n")
        fout.write("\n")

        fout.write("K_POINTS automatic\n")
        fout.write("1 1 1 0 0 0\n")
        fout.write("\n")
    xyz.to_pwscf(inp_name)


# run the simulation
for i in range(n_test + 1):
    ecut_wfc = int(ecut_wfc_min + i * ecut_wfc_step)
    inp_name = "test-ecutwfc-%d.in" % ecut_wfc
    out_f_name = "test-ecutwfc-%d.out" % ecut_wfc
    os.system("pw.x < %s | tee %s" % (inp_name, out_f_name))


# analyse the result
for i in range(n_test + 1):
    ecut_wfc = int(ecut_wfc_min + i * ecut_wfc_step)
    out_f_name = "test-ecutwfc-%d.out" % ecut_wfc
    os.system("cat %s | grep '!    total energy' >> energy-ecutwfc.data" % out_f_name)

ecut_wfc_all = [ ecut_wfc_min + i * ecut_wfc_step for i in range(n_test + 1)]
energy_all = []
with open("energy-ecutwfc.data", 'r') as fin:
    for line in fin:
        energy_all.append(float(line.split()[4]))

import matplotlib.pyplot as plt

plt.plot(ecut_wfc_all, energy_all)
plt.show()
