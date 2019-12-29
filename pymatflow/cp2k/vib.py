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
from pymatflow.cp2k.base.vibrational_analysis import cp2k_vibrational_analysis

"""
Usage:
"""

class vib_run:
    """
    Note:
        vib_run is the calss as an agent for Vibrational Analysis running.
    """
    def __init__(self):
        self.glob = cp2k_glob()
        self.force_eval = cp2k_force_eval()
        self.vibrational_analysis = cp2k_vibrational_analysis()

        self.glob.basic_setting(run_type="VIBRATIONAL_ANALYSIS")
        self.force_eval.basic_setting()
        # calculation of IR through vib need print dipole moments
        # throught DFT/PRINT/MOMENTS
        #self.force_eval.dft.printout.print_moments()
        self.force_eval.dft.printout.moments.status == True
    
    def get_xyz(self, xyzfile):
        """
        xyz_f:
            a modified xyz formatted file(the second line specifies the cell of the 
            system).
        """
        self.force_eval.subsys.xyz.get_xyz(xyzfile)
    
    def set_params(self, force_eval={}, vibrational={}):
        """
        force_eval:
            allowing control of FORCE_EVAL/... parameters by user
        vibrational:
            allowing control of VIBRATIONAL_ANALYSIS/... parameters by user
        """

        self.force_eval.set_params(force_eval)
        self.vibrational_analysis.set_params(vibrational)


    def vib(self, directory="tmp-cp2k-vib", inpname="vib.inp", output="vib.out", 
            mpi="", runopt="gen"):
        """
        directory:
            wheere the calculation will happen
        inpname:
            input filename for the cp2k
        output:
            output filename for the cp2k
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.vibrational_analysis.to_input(fout)
        
            # gen server job comit file
            self.gen_yh(cmd="cp2k.popt", inpname=inpname, output=output)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s cp2k.psmp -in %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")
   # 

    def gen_yh(self,inpname, output, directory="tmp-cp2k-vib", cmd="cp2k.psmp"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s -in %s | tee %s\n" % (cmd, inpname, output))
