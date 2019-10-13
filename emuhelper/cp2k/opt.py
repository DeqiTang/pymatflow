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

        cutoff = 60
        rel_cutoff = 30

        
    def geo_opt(self, directory="tmp-cp2k-geo-opt", inpname="geo-opt.inp", output="geo-opt.out", mpi="", runopt="gen"):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))

            self.set_geo_opt()
            self.glob.to_input(os.path.join(directory, inpname))
            self.force_eval.to_input(os.path.join(directory, inpname))
            self.motion.to_input(os.path.join(directory, inpname))
        
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s cp2k.psmp -in %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")
    
    def cell_opt(self, directory="tmp-cp2k-cell-opt", inpname="cell-opt.inp", output="cell-opt.out", mpi="", runopt="gen"):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))
            
            self.set_cell_opt()
            self.glob.to_input(os.path.join(directory, inpname))
            self.force_eval.to_input(os.path.join(directory, inpname))
            self.motion.to_input(os.path.join(directory, inpname))
        
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s cp2k.psmp -in %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def analysis(self, directory="tmp-cp2k-geo-opt", output="geo-opt.out"):
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
    
    def set_geo_opt(self):
        self.run_type = "GEO_OPT"
        self.glob.params["RUN_TYPE"] = "GEO_OPT"
        self.motion.set_type("GEO_OPT")
    
    def set_cell_opt(self):
        # Note:
        #   if you are doing CELL_OPT run, you must also enable
        #   "STRESS_TENSOR" in FORCE_EVAL%STRESS_TENSOR
        self.run_type = "CELL_OPT" 
        self.glob.params["RUN_TYPE"] = "CELL_OPT"
        self.motion.set_type("CELL_OPT")
        if self.force_eval.params["STRESS_TENSOR"] is None:
            self.force_eval.params["STRESS_TENSOR"] = "ANALYTICAL" # NUMERICAL
