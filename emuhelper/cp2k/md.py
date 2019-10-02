#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from emuhelper.cp2k.base.glob import cp2k_glob
from emuhelper.cp2k.base.force_eval import cp2k_force_eval
from emuhelper.cp2k.base.motion import cp2k_motion

"""
Usage:
"""

class md_run:
    """
    """
    def __init__(self, xyz_f):
        self.glob = cp2k_glob()
        self.force_eval = cp2k_force_eval(xyz_f)
        self.motion = cp2k_motion()

        cutoff = 60
        rel_cutoff = 30
        self.glob.params["RUN_TYPE"] = "MD"

        self.motion.set_type("MD")
        self.motion.md.params["STEPS"] = 20

    def gen_input(self, directory="tmp-md-cp2k", inpname="molecular-dynamics.inp"):
        """
        directory: a place for all the generated files
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))

        self.glob.to_input(os.path.join(directory, inpname))
        self.force_eval.to_input(os.path.join(directory, inpname))
        self.motion.to_input(os.path.join(directory, inpname))
    
    def run(self, directory="tmp-md-cp2k", inpname="molecular-dynamics.inp", output="molecular-dynamics.out"):
        """
        directory: a place for all the generated files
        """
        os.chdir(directory)
        os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
        os.chdir("../")

    def analysis(self, directory="tmp-md-cp2k", output="molecular-dynamics.out"):
        # analyse the result
        os.chdir(directory)
        os.system("cat %s | grep 'ENERGY| Total FORCE_EVAL' > energy-per-ion-step.data" % (output))
        energies = []
        with open("energy-per-ion-step.data", 'r') as fin:
            for line in fin:
                energies.append(float(line.split()[8]))

        steps = [i for i in range(len(energies))]
        plt.plot(steps, energies)
        plt.show()

        os.chdir("../")

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

