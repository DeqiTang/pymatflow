#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import matplotlib.pyplot as plt

from emuhelper.abinit.base.electrons import abinit_electrons
from emuhelper.abinit.base.ions import abinit_ions
from emuhelper.abinit.base.system import abinit_system


class md_run:
    """
    """
    def __init__(self, xyz_f):
        self.system = abinit_system(xyz_f)
        self.electrons = abinit_electrons()
        self.ions = abinit_ions()
        
        self.electrons.params["ecut"] = 50
        self.electrons.params["kptopt"] = 1
        self.electrons.params["ngkpt"] = "1 1 1"
        self.electrons.params["occopt"] = 3  # fermi dirac smearing of occupation
        self.electrons.params["nstep"] = 100
        self.electrons.params["diemac"] = 2.0
        self.ions.params["ionmov"] = 8 #11
        self.ions.params["dtion"] = 100
        self.ions.params["ntime"] = 1000
        self.ions.params["nctime"] = 1
        self.ions.params["mdtemp(1)"] = 300 
        self.ions.params["mdtemp(2)"] = 300
        self.ions.params["tolmxf"] = 5.0e-4 # Ha/Bohr
        self.ions.params["toldfe"] = 1.0e-6
        
    def gen_input(self, directory="tmp-abinit-md", inpname="molecular-dynamics.in"):
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        os.system("cp *.psp8 %s/" % directory) 

        with open(os.path.join(directory, inpname), 'w') as fout:
            self.electrons.to_in(fout)
            self.ions.to_in(fout)
            self.system.to_in(fout)

        with open(os.path.join(directory, inpname.split(".")[0]+".files"), 'w') as fout:
            fout.write("%s\n" % inpname)
            fout.write("%s.out\n" % inpname.split(".")[0])
            fout.write("%si\n" % inpname.split(".")[0])
            fout.write("%so\n" % inpname.split(".")[0])
            fout.write("temp\n")
            for element in self.system.xyz.specie_labels:
                fout.write("%s\n" % (element + ".psp8"))
    def run(self, directory="tmp-abinit-md", inpname="molecular-dynamics.in"):
        os.chdir(directory)
        os.system("abinit < %s" % inpname.split(".")[0]+".files")
        os.chdir("../")
    
    def analysis(self, directory="tmp-abinit-md", inpname="molecular-dynamics.in"):
        pass
