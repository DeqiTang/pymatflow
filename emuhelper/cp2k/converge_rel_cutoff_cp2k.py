#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from emuhelper.cp2k.base.xyz import cp2k_xyz

"""
Usage:
    python converge_rel_cutoff.py xxx.xyz rel_cutoff_min rel_cutoff_max rel_cutoff_step cutoff
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.
"""




# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

rel_cutoff_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
rel_cutoff_max = int(sys.argv[3])
rel_cutoff_step = int(sys.argv[4])
cutoff = int(sys.argv[5])

xyz = cp2k_xyz(sys.argv[1])

base_project_name = "test"

if os.path.exists("./tmp-rel-cutoff"):
    shutil.rmtree("./tmp-rel-cutoff")
os.mkdir("./tmp-rel-cutoff")
os.chdir("./tmp-rel-cutoff")
shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")

n_test = int((rel_cutoff_max - rel_cutoff_min) / rel_cutoff_step)
for i in range(n_test + 1):
    rel_cutoff = int(rel_cutoff_min + i * rel_cutoff_step)
    inp_name = "test-rel-cutoff-%d.inp" % rel_cutoff
    with open(inp_name, 'w') as fout:
        fout.write("&GLOBAL\n")
        fout.write("\tPROJECT\t%s\n" % (base_project_name + str(rel_cutoff)))
        fout.write("\tRUN_TYPE ENERGY_FORCE\n")
        fout.write("\tPRINT_LEVEL LOW\n")
        fout.write("&END GLOBAL\n")
        fout.write("\n")

        fout.write("&FORCE_EVAL\n")
        fout.write("\tMETHOD Quickstep\n")
    # subsys
    xyz.to_subsys(inp_name)
    # end subsys
    with open(inp_name, 'a') as fout:
        # dft
        fout.write("\t&DFT\n")
        fout.write("\t\tBASIS_SET_FILE_NAME BASIS_SET\n")
        fout.write("\t\tPOTENTIAL_FILE_NAME GTH_POTENTIALS\n")
        fout.write("\t\t&QS\n")
        fout.write("\t\t\tEPS_DEFAULT 1.0E-10\n")
        fout.write("\t\t&END QS\n")
        fout.write("\t\t&MGRID\n")
        fout.write("\t\t\tNGRIDS 4\n")
        fout.write("\t\t\tCUTOFF %d\n" % cutoff)
        fout.write("\t\t\tREL_CUTOFF %d\n" % rel_cutoff)
        fout.write("\t\t&END MGRID\n")
        fout.write("\t\t&XC\n")
        fout.write("\t\t\t&XC_FUNCTIONAL PADE\n")
        fout.write("\t\t\t&END XC_FUNCTIONAL\n")
        fout.write("\t\t&END XC\n")
        fout.write("\t\t&SCF\n")
        fout.write("\t\t\tSCF_GUESS ATOMIC\n")
        fout.write("\t\t\tEPS_SCF 1.0E-06\n")
        fout.write("\t\t\tMAX_SCF 300\n")
        fout.write("\t\t\t&DIAGONALIZATION ON\n")
        fout.write("\t\t\t\tALGORITHM STANDARD\n")
        fout.write("\t\t\t&END DIAGONALIZATION\n")
        fout.write("\t\t\t&MIXING T\n")
        fout.write("\t\t\t\tMETHOD BROYDEN_MIXING\n")
        fout.write("\t\t\t\tALPHA 0.4\n")
        fout.write("\t\t\t\tNBROYDEN 8\n")
        fout.write("\t\t\t&END MIXING\n")
        fout.write("\t\t&END SCF\n")
        fout.write("\t&END DFT\n")
        # end dft
        fout.write("&END FORCE_EVAL\n")

# run the simulation
for i in range(n_test + 1):
    rel_cutoff = int(rel_cutoff_min + i * rel_cutoff_step)
    inp_name = "test-rel-cutoff-%d.inp" % rel_cutoff
    out_f_name = "test-rel-cutoff-%d.out" % rel_cutoff
    os.system("cp2k.psmp -in %s | tee %s" % (inp_name, out_f_name))


# analyse the result
for i in range(n_test + 1):
    rel_cutoff = int(rel_cutoff_min + i * rel_cutoff_step)
    out_f_name = "test-rel-cutoff-%d.out" % rel_cutoff
    os.system("cat %s | grep 'Total energy:' >> energy-rel-cutoff.data" % out_f_name)

rel_cut = [rel_cutoff_min + i * rel_cutoff_step for i in range(n_test + 1)]
energy = []
with open("energy-rel-cutoff.data", 'r') as fin:
    for line in fin:
        energy.append(float(line.split()[2]))

import matplotlib.pyplot as plt

plt.plot(rel_cut, energy)
plt.show()
