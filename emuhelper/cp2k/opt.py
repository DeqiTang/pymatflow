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

class opt_run:
    """
    """
    def __init__(self, xyz_f):
        self.glob = cp2k_glob()
        self.force_eval = cp2k_force_eval(xyz_f)
        self.motion = cp2k_motion()
        
        self.run_type = "GEO_OPT" # CELL_OPT
        self.set_run_type("GEO_OPT")

        cutoff = 60
        rel_cutoff = 30

        
    def gen_input(self, directory="tmp-cp2k-opt", inpname="geometric-optimization.inp"):
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
    
    def run(self, directory="tmp-cp2k-opt", inpname="geometric-optimization.inp", output="geometric-optimization.out"):
        """
        directory: a place for all the generated files
        """
        os.chdir(directory)
        os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
        os.chdir("../")    

    def analysis(self, directory="tmp-cp2k-opt", output="geometric-optimization.out"):
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
    
    def set_run_type(self, run_type="GEO_OPT"):
        # run_type can only be one of "GEO_OPT" and "CELL_OPT"
        # Note:
        #   if you are doing CELL_OPT run, you must also enable
        #   "STRESS_TENSOR" in FORCE_EVAL%STRESS_TENSOR
        if run_type != "GEO_OPT" and run_type != "CELL_OPT":
            print("==========================================\n")
            print("           WARNING    !!!!!!!!!\n")
            print("==========================================\n")
            print("cp2k.opt.opt_run can only conduct 'GEO_OPT'\n")
            print("and 'CELL_OPT'\n")
            sys.exit(1)
        self.run_type = run_type 
        self.glob.params["RUN_TYPE"] = run_type
        self.motion.set_type(run_type)
        if run_type == "CELL_OPT" and self.force_eval.params["STRESS_TENSOR"] is None:
            self.force_eval.params["STRESS_TENSOR"] = "ANALYTICAL" # NUMERICAL
