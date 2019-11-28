#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil
import matplotlib.pyplot as plt

from pymatflow.abinit.base.dfpt import abinit_dfpt
from pymatflow.abinit.base.electrons import abinit_electrons
from pymatflow.abinit.base.system import abinit_system
from pymatflow.abinit.base.properties import abinit_properties

class dfpt_run():
    """
    procedure for DFPT calculation:
        1) ground-state calculation
        2) DFPT calculation
        3) post-processing
    Reference:
        https://docs.abinit.org/guide/respfn/
        https://docs.abinit.org/tutorial/elastic/
    """
    def __init__(self, xyz_f):
        self.system = abinit_system(xyz_f)
        self.electrons = abinit_electrons()
        self.properties = abinit_properties()
        self.dfpt = abinit_dfpt()

        self.electrons.basic_setting()
        self.dfpt.basic_setting()
        
        self.electrons.kpoints.params["kptopt"] = 2
        
    def do_dfpt(self, directory="tmp-abinit-static", inpname="dfpt.in", mpi="", runopt="gen",
            electrons={}, kpoints={}, properties=[]):
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("dfpt calculation:\n")
            print("  directory of previous static calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":

            self.electrons.set_scf_nscf("scf")
            self.electrons.set_params(electrons)
            self.electrons.kpoints.set_params(kpoints)
            self.properties.get_option(option=properties)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.electrons.to_in(fout)
                self.dfpt.to_in(fout)
                self.properties.to_in(fout)
                self.system.to_in(fout)

            with open(os.path.join(directory, inpname.split(".")[0]+".files"), 'w') as fout:
                fout.write("%s\n" % inpname)
                fout.write("%s.out\n" % inpname.split(".")[0])
                fout.write("%s-output\n" % "static-scf")
                fout.write("%s-output\n" % inpname.split(".")[0])
                fout.write("temp\n")
                for element in self.system.xyz.specie_labels:
                    fout.write("%s\n" % (element + ".psp8"))
                    #fout.write("%s\n" % (element + ".GGA_PBE-JTH.xml"))
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.chdir("../")
    # 
