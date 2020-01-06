#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import matplotlib.pyplot as plt

from pymatflow.abinit.abinit import abinit

class md_run(abinit):
    """
    """
    def __init__(self):
        super().__init__()
        
        self.electrons.basic_setting()

        self.ions.basic_setting(mode="md")

        self.guard.set_queen(queen="md", electrons=self.electrons, ions=self.ions, system=self.system)
        

    def md(self, directory="tmp-abinit-md", inpname="molecular-dynamics.in", mpi="", runopt="gen"):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory) 
            os.system("cp %s %s/" % (self.system.xyz.file, directory))
            
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
