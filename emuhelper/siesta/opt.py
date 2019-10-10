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
        
        self.electrons.xc["functional"] = "GGA"
        self.electrons.xc["authors"] = "PBE"
        self.electrons.dm["Tolerance"] = "1.d-6"
        self.electrons.dm["MixingWight"] = 0.1
        self.electrons.dm["NumberPulay"] = 5
        self.electrons.dm["AllowExtrapolation"] = "true"
        self.electrons.dm["UseSaveDM"] = "true"
        self.electrons.params["SolutionMethod"] = "diagon"
        self.electrons.params["MeshCutoff"] = 300 #100
        
        self.ions.md["TypeOfRun"] = "CG"   # CG, Broyden, 
        self.ions.md["VariableCell"] = "false"
        self.ions.md["ConstantVolume"] = "true"
        self.ions.md["MaxForceTol"] = 0.001 # eV/Ang
        self.ions.md["MaxStressTol"] = 0.01 # GPa
        self.ions.md["Steps"] = 60
        self.ions.md["MaxDispl"] = 0.2 # Bohr
        self.ions.md["PreconditionVariableCell"] = 5 # Ang

        self.ions.params["WriteCoorXmol"] = "true"
        self.ions.params["WriteMDXmol"] = "true"

    def gen_input(self, directory="tmp-siesta-opt", inpname="geometric-optimization.fdf"):
        
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        
        for element in self.system.xyz.specie_labels:
            shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))

                
        with open(os.path.join(directory, inpname), 'w') as fout:
            self.system.to_fdf(fout)
            self.electrons.to_fdf(fout)
            self.ions.to_fdf(fout)
    
    def run(self, directory="tmp-siesta-opt", inpname="geometric-optimization.fdf", output="geometric-optimization.out"):
        # run the simulation
        os.chdir(directory)
        os.system("siesta < %s | tee %s" % (inpname, output))
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
