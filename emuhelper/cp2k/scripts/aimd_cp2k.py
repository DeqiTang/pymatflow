#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from emuhelper.cp2k.base.glob import cp2k_glob
from emuhelper.cp2k.base.force_eval import cp2k_force_eval
from emuhelper.cp2k.base.motion import cp2k_motion

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

if __name__ == '__main__':
    glob = cp2k_glob()
    force_eval = cp2k_force_eval(sys.argv[1])
    motion = cp2k_motion()

    if os.path.exists("./tmp-aimd"):
        shutil.rmtree("./tmp-aimd")
    os.mkdir("./tmp-aimd")
    os.chdir("./tmp-aimd")
    shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])

    inp_name = "aimd.inp"
    glob.params["RUN_TYPE"] = "MD"
    glob.to_global(inp_name)
    force_eval.to_force_eval(inp_name)

    motion.set_type("MD")
    motion.md.params["STEPS"] = 20
    motion.to_motion(inp_name)
    
    #fout.write("\t&MD\n")
    #fout.write("\t\tENSEMBLE NVT\n")
    #fout.write("\t\tSTEPS 100\n")
    #fout.write("\t\tTIMESTEP 0.5\n")
    #fout.write("\t\t&THERMOSTAT\n")
    #fout.write("\t\t\tTYPE NOSE\n")
    #fout.write("\t\t\t&NOSE\n")
    #fout.write("\t\t\t\tTIMECON 10.0\n")
    #fout.write("\t\t\t&END NOSE\n")
    #fout.write("\t\t&END THERMOSTAT\n")
    #fout.write("\t\tTEMPERATURE 300.0\n")
    #fout.write("\t&END MD\n")
    #fout.write("\t&PRINT\n")
    #fout.write("\t\t&RESTART\n")
    #fout.write("\t\t\t&EACH\n")
    #fout.write("\t\t\t\tMD 0\n")
    #fout.write("\t\t\t&END EACH\n")
    #fout.write("\t\t&END RESTART\n")
    #fout.write("\t&END PRINT\n")
    #fout.write("&END MOTION\n")

    # run the simulation
    out_f_name = "aimd.out"
    os.system("cp2k.psmp -in %s | tee %s" % (inp_name, out_f_name))

    # analyse the result

    import matplotlib.pyplot as plt

