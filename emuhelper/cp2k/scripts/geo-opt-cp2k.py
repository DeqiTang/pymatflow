#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import os
import shutil

import matplotlib.pyplot as plt

from emuhelper.cp2k.base.glob import cp2k_glob
from emuhelper.cp2k.base.force_eval import cp2k_force_eval
from emuhelper.cp2k.base.motion import cp2k_motion

"""
Usage:
    python geo_opt_cp2k.py xxx.xyz 
    xxx.xyz is the input structure file

    make sure the xyz structure file 
    is in the directory.
"""

def analysis_result(out_f_name):
    # analyse the result
    os.system("cat %s | grep 'Total Energy' > energy-per-geo-step.data" % (out_f_name))
    energies = []
    with open("energy-per-geo-step.data", 'r') as fin:
        for line in fin:
            energies.append(float(line.split()[3]))

    steps = [i for i in range(len(energies))]
    plt.plot(steps, energies)
    plt.show()

if __name__ == "__main__":
    inp_name = "geo-opt.inp"
    out_f_name = "geo-opt.out"
    glob = cp2k_glob()
    force_eval = cp2k_force_eval(sys.argv[1])
    motion = cp2k_motion()

    if os.path.exists("tmp-geo-opt"):
        shutil.rmtree("tmp-geo-opt")
    os.mkdir("tmp-geo-opt")
    os.chdir("tmp-geo-opt")
    shutil.copyfile("../%s" % sys.argv[1], "./%s" % sys.argv[1])

    glob.params["RUN_TYPE"] = "GEO_OPT"
    glob.to_global(inp_name)
    force_eval.to_force_eval(inp_name)
    motion.set_type("GEO_OPT")
    motion.to_motion(inp_name)
    os.system("cp2k.psmp -in %s | tee log.txt" % inp_name)
    analysis_result(out_f_name)



