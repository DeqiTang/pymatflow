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

"""
Usage:
"""

class static_run:
    """
    """
    def __init__(self, xyz_f):
        self.glob = cp2k_glob()
        self.force_eval = cp2k_force_eval(xyz_f)

        self.glob.params["RUN_TYPE"] = "ENERGY_FORCE"
        self.force_eval.dft.mgrid.params["CUTOFF"] = 100
        self.force_eval.dft.mgrid.params["REL_CUTOFF"]= 60


    def gen_input(self, directory="tmp-cp2k-static", inpname="static.inp"):
        """
        directory: a place for all the generated files
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))

        self.glob.to_input(os.path.join(directory, inpname))
        self.force_eval.to_input(os.path.join(directory, inpname))
    
    def run(self, directory="tmp-cp2k-static", inpname="static.inp", output="static.out"):
        """
        directory: a place for all the generated files
        """
        os.chdir(directory)
        os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
        os.chdir("../")    

    def analysis(self, directory="tmp-cp2k-static", output="static.out"):
        # analyse the result
        os.chdir(directory)

        os.chdir("../")
    
    def converge_cutoff(self, emin, emax, step, rel_cutoff, directory="tmp-cp2k-cutoff"):
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))
        
        n_test = int((emax - emin) / step)
        for i in range(n_test + 1):
            cutoff = int(emin + i * step)
            inpname = "cutoff-%d" % cutoff
            self.force_eval.dft.mgrid.params["CUTOFF"] = cutoff
            self.force_eval.dft.mgrid.params["REL_CUTOFF"] = rel_cutoff
            self.glob.to_input(os.path.join(directory, inpname))
            self.force_eval.to_input(os.path.join(directory, inpname))
        # run
        os.chdir(directory)
        for i in range(n_test + 1):
            cutoff = int(emin + i * step)
            inpname = "cutoff-%d" % cutoff
            output = "cutoff-%d.out" % cutoff
            os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
        os.chdir("../")
        # analyse the result
        os.chdir(directory)
        for i in range(n_test + 1):
            cutoff = int(emin + i * step)
            out_f_name = "cutoff-%d.out" % cutoff
            os.system("cat %s | grep 'Total energy:' >> energy-cutoff.data" % out_f_name)
        cutoffs = [emin + i * step for i in range(n_test + 1)]
        energy = []
        with open("energy-cutoff.data", 'r') as fin:
            for line in fin:
                energy.append(float(line.split()[2]))
        plt.plot(cutoffs, energy)
        plt.show()
        os.chdir("../")

    def converge_rel_cutoff(self, emin, emax, step, cutoff, directory="tmp-cp2k-rel-cutoff"):
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))
        
        n_test = int((emax - emin) / step)
        for i in range(n_test + 1):
            rel_cutoff = int(emin + i * step)
            inpname = "rel-cutoff-%d" % rel_cutoff
            self.force_eval.dft.mgrid.params["CUTOFF"] = cutoff
            self.force_eval.dft.mgrid.params["REL_CUTOFF"] = rel_cutoff
            self.glob.to_input(os.path.join(directory, inpname))
            self.force_eval.to_input(os.path.join(directory, inpname))
        # run
        os.chdir(directory)
        for i in range(n_test + 1):
            rel_cutoff = int(emin + i * step)
            inpname = "rel-cutoff-%d" % rel_cutoff
            output = "rel-cutoff-%d.out" % rel_cutoff
            os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
        os.chdir("../")
        # analyse the result
        os.chdir(directory)
        for i in range(n_test + 1):
            rel_cutoff = int(emin + i * step)
            out_f_name = "rel-cutoff-%d.out" % rel_cutoff
            os.system("cat %s | grep 'Total energy:' >> energy-rel-cutoff.data" % out_f_name)
        rel_cutoffs = [emin + i * step for i in range(n_test + 1)]
        energy = []
        with open("energy-rel-cutoff.data", 'r') as fin:
            for line in fin:
                energy.append(float(line.split()[2]))
        plt.plot(rel_cutoffs, energy)
        plt.show()
        os.chdir("../")

    def print_electron_density(self):
        self.force_eval.dft.printout.print_electron_density()
    
    def print_bands(self):
        self.force_eval.dft.printout.print_bands()

    def print_pdos(self):
        self.force_eval.dft.printout.print_pdos()
