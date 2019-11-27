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
    Note:
        lr_run is the  class as an agent for Linear Response calculation.
    """
    def __init__(self, xyz_f):
        self.glob = cp2k_glob()
        self.force_eval = cp2k_force_eval(xyz_f)
 #       self.atom = cp2k_atom()
        
        self.glob.basic_setting(run_type="LINEAR_RESPONSE")
        self.force_eval.basic_setting()


    def lr(self, directory="tmp-cp2k-lr", inpname="lr.inp", output="lr.out", 
            force_eval={}, mpi="", runopt="gen", printout_option=[]):
        """
        directory:
            where the calculation will happen
        inpname:
            inputfile name for the cp2k
        output:
            output filename for the cp2k
        force_eval:
            allowing control of FORCE_EVAL/... parameters by user
        printout_option:
            a list of integers, controlling the printout of properties, etc.
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
            self.gen_yh(cmd="cp2k.popt", directory=directory, inpname=inpname, output=output)   

        if runopt == "run" or runopt == "genrun":
           os.chdir(directory)
           os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
           os.chdir("../")    


    def printout_option(self, option=[]):
        """
        Note:
            responsible for the parseing of the printout_option like in self.scf()

        option:
            1: printout pdos
            2: printout band
            3: printout electron densities
            4: printout electron local function(ELF)
            5: printout molecular orbitals
            6: printout molecular orbital cube files
            7: printout mulliken populaltion analysis
            8: printout cubes for generation of STM images
            9: printout cube file with total density(electrons+atomic core)
           10: printout v_hartree_cube
           11: printout v_xc_cube
           12: printout xray_diffraction_spectrum
           13: request a RESP fit of charges.
           14: request a LINRES calculation
        """
        if 1 in option:
            self.force_eval.dft.printout.print_pdos()
        if 2 in option:
            self.force_eval.dft.printout.print_band(self.force_eval.subsys.xyz)
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
        if 14 in option:
            self.force_eval.properties.linres.status = True
            # XC_DERIV method not implemented for GPW!
            # so we try GAPW
            self.force_eval.dft.qs.params["METHOD"] = "RIGPW"

    def gen_yh(self,inpname, output, directory="tmp-cp2k-lr", cmd="cp2k.psmp"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s -in %s | tee %s\n" % (cmd, inpname, output))
