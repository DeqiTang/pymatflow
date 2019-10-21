#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


from emuhelper.siesta.base.system import siesta_system
from emuhelper.siesta.base.electrons import siesta_electrons
from emuhelper.siesta.base.ions import siesta_ions


class opt_run:
    """
    """
    def __init__(self, xyz_f):
        self.system = siesta_system(xyz_f)
        self.electrons = siesta_electrons()
        self.ions = siesta_ions()

        self.electrons.basic_setting()
        self.ions.basic_setting()

    def opt(self, directory="tmp-siesta-opt", inpname="geometric-optimization.fdf", output="geometric-optimization.out",
            mpi="", runopt="gen", mode=0, electrons={}, kpoints_mp=[1, 1, 1]):
        
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
        
            for element in self.system.xyz.specie_labels:
                shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))

            self.electrons.kpoints_mp = kpoints_mp
            self.electrons.set_params(electrons)
            self.set_opt_mode(mode=mode)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.system.to_fdf(fout)
                self.electrons.to_fdf(fout)
                self.ions.to_fdf(fout)
        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            os.system("%s siesta < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")


    def analysis(self, directory="tmp-siesta-opt", inpname="geometric-optimization.fdf", output="geometric-optimization.out"):
        # analyse the results
        import matplotlib.pyplot as plt

        os.chdir(directory)
        os.system("cat %s | grep 'siesta: E_KS(eV) =' > energy-per-ion-step.data" % (output))

        energies = []
        with open("energy-per-ion-step.data", 'r') as fin:
            for line in fin:
                energies.append(float(line.split()[3]))

        steps = [i for i in range(len(energies))]
        plt.plot(steps, energies)
        plt.show()
        os.chdir("../")

    def set_opt_mode(self, mode=0):
        """
        mode:
            0: do not variable cell
            1: variable cell
        """
        if mode == 0:
            self.ions.md["VariableCell"] = "false"
        elif mode == 1:
            self.ions.md["VariableCell"] = "false"
        else:
            print("==========================================\n")
            print("             WARNING !!!\n")
            print("opt mode can only be 0 or 1\n")
            print("where 0 is: MD.VariableCell = flase\n")
            print("and 1 is : MD.VariableCell = true\n")
            sys.exit(1)
