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
    python aimd_cp2k.py xxx.xyz 
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.
"""

# =======================================
#           Control Parameters
# =======================================

cutoff = 60

rel_cutoff = 30

xyz = cp2k_xyz(sys.argv[1])

base_project_name = "aimd"

potential_file = "POTENTIAL"
basis_file = "BASIS_SET"


# build folder to conduct the computing

if os.path.exists("./tmp-aimd"):
    shutil.rmtree("./tmp-aimd")
os.mkdir("./tmp-aimd")
os.chdir("./tmp-aimd")
shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])

inp_name = "aimd.inp"
with open(inp_name, 'w') as fout:
    fout.write("&GLOBAL\n")
    fout.write("\tPROJECT\t%s\n" % (base_project_name ))
    fout.write("\tRUN_TYPE MD\n")
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
    fout.write("\t\tBASIS_SET_FILE_NAME %s\n" % basis_file)
    fout.write("\t\tPOTENTIAL_FILE_NAME %s\n" % potential_file)
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
    fout.write("\t&MD\n")
    fout.write("\t\tENSEMBLE NVT\n")
    fout.write("\t\tSTEPS 100\n")
    fout.write("\t\tTIMESTEP 0.5\n")
    fout.write("\t\t&THERMOSTAT\n")
    fout.write("\t\t\tTYPE NOSE\n")
    fout.write("\t\t\t&NOSE\n")
    fout.write("\t\t\t\tTIMECON 10.0\n")
    fout.write("\t\t\t&END NOSE\n")
    fout.write("\t\t&END THERMOSTAT\n")
    fout.write("\t\tTEMPERATURE 300.0\n")
    fout.write("\t&END MD\n")
    fout.write("\t&PRINT\n")
    fout.write("\t\t&RESTART\n")
    fout.write("\t\t\t&EACH\n")
    fout.write("\t\t\t\tMD 0\n")
    fout.write("\t\t\t&END EACH\n")
    fout.write("\t\t&END RESTART\n")
    fout.write("\t&END PRINT\n")
    fout.write("&END MOTION\n")
    fout.write("\n")

# run the simulation
out_f_name = "aimd.out"
os.system("cp2k.psmp -in %s | tee %s" % (inp_name, out_f_name))


# analyse the result

import matplotlib.pyplot as plt

