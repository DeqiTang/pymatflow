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
#from emuhelper.cp2k.base.atom import cp2k_atom

"""
Usage:
"""

class lr_run:
    """
    """
    def __init__(self, xyz_f):
        self.glob = cp2k_glob()
        self.force_eval = cp2k_force_eval(xyz_f)
 #       self.atom = cp2k_atom()
        
        self.glob.basic_setting(run_type="LINEAR_RESPONSE")
        self.force_eval.basic_setting()


    def lr(self, directory="tmp-cp2k-lr", inpname="lr.inp", output="lr.out", 
            force_eval={}, mpi="", runopt="gen", printout_option=0):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))
            # using force_eval
            self.force_eval.set_params(force_eval)
            #self.atom.set_params(atom)
            self.printout_option(printout_option)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                #self.atom.to_input(fout)
 
            # gen server job comit file
            self.gen_yh(cmd="cp2k.popt", inpname=inpname, output=output)   

        if runopt == "run" or runopt == "genrun":
           os.chdir(directory)
           os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
           os.chdir("../")    


    def printout_option(self, option=[]):
        """
        option:
            0: do not printout properties
            1: printout pdos
            2: printout bands
            3: printout electron densities
        """
        if 1 in option:
            self.force_eval.dft.printout.print_pdos()
        if 2 in option:
            self.force_eval.dft.printout.print_bands()
        if 3 in option:
            self.force_eval.dft.printout.print_electron_density()
        if 4 in option:
            self.force_eval.dft.printout.elf_cube = True
        if 5 in option:
            self.force_eval.dft.printout.mo = True
        if 6 in option:
            self.force_eval.dft.printout.mo_cubes = True
        if 7 in option:
            self.force_eval.dft.printout.mulliken = True
        if 8 in option:
            self.force_eval.dft.printout.stm = True
        if 9 in option:
            self.force_eval.dft.printout.tot_density_cube = True
        if 10 in option:
            self.force_eval.dft.printout.v_hartree_cube = True
        if 11 in option:
            self.force_eval.dft.printout.v_xc_cube = True
        if 12 in option:
            self.force_eval.dft.printout.xray_diffraction_spectrum = True
        if 13 in option:
            self.force_eval.properties.resp.status = True


    def gen_yh(self,inpname, output, directory="tmp-cp2k-static", cmd="cp2k.psmp"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s -in %s | tee %s\n" % (cmd, inpname, output))
