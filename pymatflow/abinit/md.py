#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import matplotlib.pyplot as plt

from pymatflow.abinit.base.electrons import abinit_electrons
from pymatflow.abinit.base.ions import abinit_ions
from pymatflow.abinit.base.system import abinit_system

from pymatflow.abinit.base.guard import abinit_guard

class md_run:
    """
    """
    def __init__(self, xyz_f):
        self.system = abinit_system(xyz_f)
        self.electrons = abinit_electrons()
        self.ions = abinit_ions()
        
        self.electrons.basic_setting()

        self.ions.basic_setting(mode="md")

        self.guard = abinit_guard(queen="md", electrons=self.electrons, ions=self.ions, system=self.system)
        
    def md(self, directory="tmp-abinit-md", inpname="molecular-dynamics.in", mpi="", runopt="gen",
            electrons={}, ions={}, kpoints={}):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory) 
            os.system("cp %s %s/" % (self.system.xyz.file, directory))
            
            self.electrons.set_params(electrons)
            self.electrons.kpoints.set_params(kpoints)
            self.ions.set_params(ions)
            #
            self.electrons.params["tolvrs"] = 1.0e-10
            self.electrons.params["tolwrf"] = None
            self.electrons.params["toldff"] = None
            self.electrons.params["tolrff"] = None
            self.electrons.params["toldfe"] = None #1.0e-6
            #
            self.guard.check_all()
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
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.chdir("../")
    
    def analysis(self, directory="tmp-abinit-md", inpname="molecular-dynamics.in"):
        pass
