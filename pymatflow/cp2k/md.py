#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.cp2k.base.glob import cp2k_glob
from pymatflow.cp2k.base.force_eval import cp2k_force_eval
from pymatflow.cp2k.base.motion import cp2k_motion

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

        self.glob.basic_setting(run_type="MD")
        self.force_eval.basic_setting()

        self.motion.set_type("MD")

    def md(self, directory="tmp-cp2k-md", inpname="md.inp", output="md.out", mpi="", runopt="gen",
            force_eval={}, motion={}):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))

            self.force_eval.set_params(force_eval)
            self.motion.set_params(motion)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.motion.to_input(fout)

            # gen server job comit file
            self.gen_yh(cmd="cp2k.popt", inpname=inpname, output=output)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s cp2k.psmp -in %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")
    

    def analysis(self, directory="tmp-cp2k-md", output="md.out"):
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

    def ir_spectra(self):
        """
        if you are calculating ir spectra, you have to
        have the dipole information for the molecule 
        available in the simulated trajectory.
        that is realized by FORCE_EVAL%DFT%LOCALIZE

        Reference:
            http://www.travis-analyzer.de/
        """
        self.force_eval.dft.localize.status = True

    def gen_yh(self,inpname, output, directory="tmp-cp2k-static", cmd="cp2k.psmp"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s -in %s | tee %s\n" % (cmd, inpname, output))
