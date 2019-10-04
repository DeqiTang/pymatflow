#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import re
import pymatgen as mg

from emuhelper.qe.base.control import qe_control
from emuhelper.qe.base.system import qe_system
from emuhelper.qe.base.electrons import qe_electrons
from emuhelper.qe.base.ions import qe_ions
from emuhelper.qe.base.arts import qe_arts

"""
Usage:
    python converge-ecutwfc-pwscf.py xxx.xyz ecut_wfc_min ecut_wfc_max ecut_wfc_step
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


ctl = qe_control()
system = qe_system()
electrons = qe_electrons()
ions = qe_ions()
arts = qe_arts(sys.argv[1])
ctl.params['title'] = "Test ecutwfc"
ctl.params['prefix'] = 'test-ecutwfc'
ctl.params['calculation'] = 'scf'
ctl.params['pseudo_dir'] = './'
system.params["ibrav"] = 0
system.params["nat"] = arts.xyz.natom
system.params["ntyp"] = arts.xyz.nspecies
system.params["input_DFT"] = "PBE"

if os.path.exists("./tmp-ecutwfc"):
    shutil.rmtree("./tmp-ecutwfc")
os.mkdir("./tmp-ecutwfc")
os.chdir("./tmp-ecutwfc")
os.system("cp ../*.UPF ./")

n_test = int((ecut_wfc_max - ecut_wfc_min) / ecut_wfc_step)
for i in range(n_test + 1):
    ecut_wfc = int(ecut_wfc_min + i * ecut_wfc_step)
    ecut_rho = ecut_wfc * 4 # using default value for ecut_rho: 4 * ecutwfc
    inp_name = "test-ecutwfc-%d.in" % ecut_wfc
    ctl.params['outdir'] = './tmp-' + str(ecut_wfc)
    system.params['ecutwfc'] = ecut_wfc
    with open(inp_name, 'w') as fout:
        ctl.to_in(fout)
        system.to_in(fout)
        electrons.to_in(fout)
        arts.to_in(fout)


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
