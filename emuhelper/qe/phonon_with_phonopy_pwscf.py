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
    python phonon_with_phonopy_pwscf.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential
    file for all the elements in the system is in the directory.

Note:
    现在只支持设置ATOMIC_POSITIONS 为crystal类型
    而我喜欢用angstrom, 所以就暂且搁置, 等待以后其支持angstrom
    参考:
    https://atztogo.github.io/phonopy/qe.html
"""



# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

#ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#ecut_max = int(sys.argv[3])
#ecut_step = int(sys.argv[4])

ecut_wfc = 35
ecut_rho = ecut_wfc * 4

supercell_n = "1 1 1"

xyz = qe_xyz(sys.argv[1])

title = "Phonon with Phonopy"

base_prefix = "bfo"

if os.path.exists("./tmp-phonon-with-phonopy"):
    shutil.rmtree("./tmp-phonon-with-phonopy")
os.mkdir("./tmp-phonon-with-phonopy")
os.chdir("./tmp-phonon-with-phonopy")

os.system("cp ../*.UPF ./")

head_inp_name = "phonon_head.in"
with open(head_inp_name, 'w') as fout:
    fout.write("&control\n")
    fout.write("calculation = 'scf'\n")
    fout.write("title = '%s'\n" % title)
    fout.write("prefix = '%s'\n" % (base_prefix))
    fout.write("restart_mode = 'from_scratch'\n")
    fout.write("nstep = 300\n")
    fout.write("outdir = '%s'\n" % ("./tmp"))
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
    fout.write("1 1 1 0 0 0\n")
    fout.write("\n")
# xyz.to_pwscf(inp_name)
pos_inp_name = "pos.in"
os.system("cat %s > %s" % (head_inp_name, pos_inp_name))
xyz.to_pwscf(pos_inp_name)

# set up the Phonopy calculation
os.system("phonopy --qe -d --dim='%s' -c %s" % (supercell_n, pos_inp_name))
os.system("ls | grep 'supercell-' > pos.data")
disp_dirs = []
with open("pos.data", 'r') as fin:
    for line in fin:
        disp_dirs.append(line.split(".")[0].split("-")[1])

for disp in disp_dirs:
    os.system("cat %s supercell-%s.in > supercell-%s-full.in" % (head_inp_name, disp, disp))
    os.system("rm supercell-%s.in" % disp)
# run the dft
for disp in disp_dirs:
    os.system("pw.x < supercell-%s-full.in > supercell-%s.out" % (disp, disp))

# analyse the result
import matplotlib.pyplot as plt

os.system("phonopy --qe -f supercell-{001..%s}.out" % (disp_dirs[-1]))

# plot band structure
with open("band.conf", 'w') as fout:
    fout.write("ATOM_NAME =")
    for element in xyz.specie_labels:
        fout.write(" %s" % element)
    fout.write("\n")
    fout.write("DIM = %s\n" % supercell_n)
    fout.write("BAND = 0.5 0.5 0.5 0.0 0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.0\n")
os.system("phonopy --qe -c %s -p band.conf" % pos_inp_name)
