#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from emuhelper.siesta.base.system import siesta_system
from emuhelper.siesta.base.electrons import siesta_electrons
from emuhelper.siesta.base.properties import siesta_properties

class static_run:
    """
    """
    def __init__(self, xyz_f):
        self.system = siesta_system(xyz_f)
        self.electrons = siesta_electrons()
        self.properties = siesta_properties()
        
        self.electrons.xc["functional"] = "GGA"
        self.electrons.xc["authors"] = "PBE"
        self.electrons.dm["Tolerance"] = "1.d-6"
        self.electrons.dm["MixingWight"] = 0.1
        self.electrons.dm["NumberPulay"] = 5
        self.electrons.dm["AllowExtrapolation"] = "true"
        self.electrons.dm["UseSaveDM"] = "false"
        self.electrons.params["SolutionMethod"] = "diagon"
        self.electrons.params["MeshCutoff"] = 100
        

    def gen_input(self, directory="tmp-static-siesta", inpname="static.fdf"):
        
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        
        for element in self.system.xyz.specie_labels:
            shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))

                
        with open(os.path.join(directory, inpname), 'w') as fout:
            self.system.to_fdf(fout)
            self.electrons.to_fdf(fout)
            self.properties.to_fdf(fout)
    
    def run(self, directory="tmp-static-siesta", inpname="static.fdf", output="static.out"):
        # run the simulation
        os.chdir(directory)
        os.system("siesta < %s | tee %s" % (inpname, output))
        os.chdir("../")


    def analysis(self, directory="tmp-static-siesta", inpname="static.fdf", output="static.out"):
        # analyse the results
        
        if self.properties.option == "pdos":
            self.dos_analysis(directory)


    def dos_analysis(self, directory):
        # plot dos
        os.chdir(directory)
        energy = []
        states = []
        with open(self.system.label+".DOS", 'r') as fin:
            for line in fin:
                energy.append(float(line.split()[0]))
                states.append(float(line.split()[1]))
        plt.plot(energy, states)
        plt.show()
        os.chdir("../")
