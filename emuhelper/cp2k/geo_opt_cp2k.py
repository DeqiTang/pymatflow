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
    python geo_opt_cp2k.py xxx.xyz 
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.
"""


        

# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

#cutoff_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#cutoff_max = int(sys.argv[3])
#cutoff_step = int(sys.argv[4])
#rel_cutoff = int(sys.argv[5])
cutoff = 60
rel_cutoff = 30

xyz = cp2k_xyz(sys.argv[1])

base_project_name = "geo-opt"

if os.path.exists("./tmp-geo-opt"):
    shutil.rmtree("./tmp-geo-opt")
os.mkdir("./tmp-geo-opt")
os.chdir("./tmp-geo-opt")
shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")

inp_name = "geo-opt.inp"
with open(inp_name, 'w') as fout:
    fout.write("&GLOBAL\n")
    fout.write("\tPROJECT\t%s\n" % (base_project_name ))
    fout.write("\tRUN_TYPE GEO_OPT\n")
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
    fout.write("\t\t\tMAX_SCF 200\n")
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

    # MOTION
    fout.write("&MOTION\n")
    fout.write("\t&GEO_OPT\n")
    fout.write("\t\tTYPE MINIMIZATION\n")
    fout.write("\t\tMAX_DR 1.0E-03\n")
    fout.write("\t\tMAX_FORCE 1.0E-03\n")
    fout.write("\t\tRMS_DR 1.0E-03\n")
    fout.write("\t\tRMS_FORCE 1.0E-03\n")
    fout.write("\t\tMAX_ITER 200\n")
    fout.write("\t\tOPTIMIZER CG\n")
    fout.write("\t\t&CG\n")
    fout.write("\t\t\tMAX_STEEP_STEPS 0\n")
    fout.write("\t\t\tRESTART_LIMIT 9.0E-01\n")
    fout.write("\t\t&END CG\n")
    fout.write("\t&END GEO_OPT\n")
    fout.write("\t&CONSTRAINT\n")
    fout.write("\t\t&FIXED_ATOMS\n")
    fout.write("\t\t\tCOMPONENTS_TO_FIX XYZ\n")
    fout.write("\t\t\tLIST\n")
    fout.write("\t\t&END FIXED_ATOMS\n")
    fout.write("\t&END CONSTRAINT\n")
    fout.write("&END MOTION\n")
    fout.write("\n")

# run the simulation
out_f_name = "geo-opt.out"
os.system("cp2k.psmp -in %s | tee %s" % (inp_name, out_f_name))


# analyse the result

import matplotlib.pyplot as plt

os.system("cat %s | grep 'Total Energy' > energy-per-geo-step.data" % (out_f_name))

energies = []
with open("energy-per-geo-step.data", 'r') as fin:
    for line in fin:
        energies.append(float(line.split()[3]))

steps = [i for i in range(len(energies))]
plt.plot(steps, energies)
plt.show()
