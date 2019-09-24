#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import re
import pymatgen as mg

from base.xyz import qe_xyz
"""
Usage:
    python converge_ecutrho_pwscf.py xxx.xyz ecut_rho_min ecut_rho_max ecut_rho_step ecut_wfc
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential
    file for all the elements in the system is in the directory.
"""


# from the xyz file: sys.argv[1]

ecut_rho_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
ecut_rho_max = int(sys.argv[3])
ecut_rho_step = int(sys.argv[4])

ecut_wfc = int(sys.argv[5])

xyz = qe_xyz(sys.argv[1])

title = "Test ecutrho"
base_prefix = "knn"

if os.path.exists("./tmp-ecutrho"):
    shutil.rmtree("./tmp-ecutrho")
os.mkdir("./tmp-ecutrho")
os.chdir("./tmp-ecutrho")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")
os.system("cp ../*.UPF ./")

n_test = int((ecut_rho_max - ecut_rho_min) / ecut_rho_step)
for i in range(n_test + 1):
    ecut_rho = int(ecut_rho_min + i * ecut_rho_step)
    inp_name = "test-ecutrho-%d.in" % ecut_rho
    with open(inp_name, 'w') as fout:
        fout.write("&control\n")
        fout.write("calculation = 'scf'\n")
        fout.write("title = '%s'\n" % title)
        fout.write("prefix = '%s'\n" % (base_prefix + str(ecut_rho)))
        fout.write("restart_mode = 'from_scratch'\n")
        fout.write("nstep = 300\n")
        fout.write("outdir = '%s'\n" % ("./tmp-" + str(ecut_rho)))
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
        fout.write("ecutrho = %d\n" % ecut_rho)
        fout.write("input_DFT = 'PBE'\n")
        fout.write("occupations = 'smearing'\n")
        fout.write("degauss = 1.0d-4\n")
        fout.write("smearing = 'marzari-vanderbilt'\n")
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
        fout.write("2 2 2 0 0 0\n")
        fout.write("\n")
    xyz.to_pwscf(inp_name)


# run the simulation
for i in range(n_test + 1):
    ecut_rho = int(ecut_rho_min + i * ecut_rho_step)
    inp_name = "test-ecutrho-%d.in" % ecut_rho
    out_f_name = "test-ecutrho-%d.out" % ecut_rho
    os.system("pw.x < %s > %s" % (inp_name, out_f_name))


# analyse the result
for i in range(n_test + 1):
    ecut_rho = int(ecut_rho_min + i * ecut_rho_step)
    out_f_name = "test-ecutrho-%d.out" % ecut_rho
    os.system("cat %s | grep '!    total energy' >> energy-ecutrho.data" % out_f_name)

ecut_rho_all = [ ecut_rho_min + i * ecut_rho_step for i in range(n_test + 1)]
energy_all = []
with open("energy-ecutrho.data", 'r') as fin:
    for line in fin:
        energy_all.append(float(line.split()[4]))

import matplotlib.pyplot as plt

plt.plot(ecut_rho_all, energy_all)
plt.show()
