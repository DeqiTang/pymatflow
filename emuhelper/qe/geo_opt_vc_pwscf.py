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
    python geo_opt_vc_pwscf.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential
    file for all the elements in the system is in the directory.
"""


# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

#ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#ecut_max = int(sys.argv[3])
#ecut_step = int(sys.argv[4])

ecut_wfc = 40 # in Ry: 1 Ry = 13.6 ev
ecut_rho = ecut_wfc * 4 # default is ecutrho = 4 * ecutwfc

xyz = qe_xyz(sys.argv[1])

title = "Geometric Optimization"

base_prefix = "geo"

if os.path.exists("./tmp-geo-opt-vc"):
    shutil.rmtree("./tmp-geo-opt-vc")
os.mkdir("./tmp-geo-opt-vc")
os.chdir("./tmp-geo-opt-vc")

os.system("cp ../*.UPF ./")

inp_name = "geo-opt-vc.in"
with open(inp_name, 'w') as fout:
    fout.write("&control\n")
    fout.write("calculation = 'vc-relax'\n")
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
    fout.write("conv_thr = 1.0d-6\n")
    fout.write("mixing_mode = 'plain'\n")
    fout.write("mixing_beta = 0.7d0\n")
    fout.write("mixing_ndim = 8") # default: 8 这个参数对于内存的使用有较大影响
    fout.write("diagonalization = 'david'\n")
    fout.write("scf_must_converge = .true.\n")
    fout.write("/\n")
    fout.write("\n")

    fout.write("&ions\n")
    fout.write("ion_dynamics = 'bfgs'\n")
    fout.write("ion_temperature ='not_controlled'\n")
    fout.write("tempw = 300.D0\n")
    fout.write("/\n")
    fout.write("\n")

    fout.write("&cell\n")
    fout.write("cell_dynamics = 'bfgs'\n")
    fout.write("press = 0.0\n")
    fout.write("/\n")
    fout.write("\n")

    fout.write("K_POINTS automatic\n")
    fout.write("1 1 1 0 0 0\n")
    fout.write("\n")
xyz.to_pwscf(inp_name)


# run the simulation
inp_name = "geo-opt-vc.in"
out_f_name = "geo-opt-vc.out"
os.system("pw.x < %s > %s" % (inp_name, out_f_name))


# analyse the result

import matplotlib.pyplot as plt
# 要注意pwscf的scf进行的次数与离子步的次数没有必然联系
# 其它程序可能scf循环的次数比离子步大一
# 但是pwscf中可能不止大一, 有些离子步进行了多次scf循环计算
os.system("cat %s | grep '!    total energy' > geo.data" % out_f_name)

energies = []
with open("geo.data", 'r') as fin:
    for line in fin:
        energies.append(float(line.split()[4]))
steps = np.arange(len(energies))
plt.plot(steps, energies)
plt.show()
