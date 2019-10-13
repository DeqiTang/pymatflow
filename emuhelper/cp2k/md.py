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

    def md(self, directory="tmp-cp2k-md", inpname="md.inp", output="md.out", mpi="", runopt="gen"):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))

            self.glob.to_input(os.path.join(directory, inpname))
            self.force_eval.to_input(os.path.join(directory, inpname))
            self.motion.to_input(os.path.join(directory, inpname))
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s cp2k.psmp < %s | tee %s" % (mpi, inpname, output))
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
