#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


from pymatflow.siesta.base.system import siesta_system
from pymatflow.siesta.base.electrons import siesta_electrons
from pymatflow.siesta.base.ions import siesta_ions


class opt_run:
    """
    """
    def __init__(self):
        self.system = siesta_system()
        self.electrons = siesta_electrons()
        self.ions = siesta_ions()

        self.electrons.basic_setting()
        self.ions.basic_setting(option="opt")

    def get_xyz(self, xyzfile):
        self.system.xyz.get_xyz(xyzfile)

    def opt(self, directory="tmp-siesta-opt", inpname="geometric-optimization.fdf", output="geometric-optimization.out",
            mpi="", runopt="gen", mode=0, electrons={}, ions={}, kpoints_mp=[1, 1, 1]):
        """
        mode:
            0: do not vary the cell
            1: vary the cell
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
        
            for element in self.system.xyz.specie_labels:
                shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))

            self.electrons.kpoints_mp = kpoints_mp
            self.electrons.set_params(electrons)
            self.ions.set_params(ions)
            self.set_opt_mode(mode=mode)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.system.to_fdf(fout)
                self.electrons.to_fdf(fout)
                self.ions.to_fdf(fout)

            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="siesta")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            os.system("%s siesta < %s | tee %s" % (mpi, inpname, output))
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

    def gen_yh(self, inpname, output, directory="tmp-siesta-static", cmd="siesta"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % (cmd, inpname, output))
