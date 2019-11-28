#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import matplotlib.pyplot as plt

from pymatflow.abinit.base.electrons import abinit_electrons
from pymatflow.abinit.base.system import abinit_system
#from pymatflow.abinit.base.ions import abinit_ions
#from emuhelper.abinit.base.properties import abinit_properties
from pymatflow.abinit.base.guard import abinit_guard

class neb_run:
    """
    """
    def __init__(self, images):
        """
        images:
            ["first.xyz", "last.xyz"]
        """
        self.electrons = abinit_electrons()
        self.system = []
        for image in images:
            self.system.append(abinit_system(image))

        self.electrons.basic_setting()
        
        self.guard = abinit_guard(queen="neb", electrons=self.electrons, system=self.system)

        self.params = {
                "imgmov": 5,
                "nimage": 7,
                "ntimimage": 30,
                "mep_solver": None,
                "tolimg": 5.0e-05,
                "imgwfstor": None,
                "mep_mxstep": None,
                "neb_algo": None,
                "neb_string": None,
                "npimage": None,
                "string_algo": None,
                "cineb_start": 7,
                "dynimage": "0 5*1 0",
                "fxcartfactor": 1.0,
                }
        
    def neb(self, directory="tmp-abinit-neb", inpname="neb.in", mpi="", runopt="gen",
            electrons={}, kpoints={}):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
            for image in self.system:
                os.system("cp %s %s/" % (image.xyz.file, directory))

            self.electrons.set_scf_nscf("scf")
            self.electrons.set_params(electrons)
            self.electrons.kpoints.set_params(kpoints)
            #
            self.guard.check_all()
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.electrons.to_in(fout)
                self.system_to_in(fout)
                fout.write("\n")
                fout.write("# ===============================\n")
                fout.write("# neb related setting\n")
                fout.write("# ===============================\n")
                for item in self.params:
                    if self.params[item] is not None:
                        fout.write("%s %s\n" % (item, str(self.params[item])))
                
            with open(os.path.join(directory, inpname.split(".")[0]+".files"), 'w') as fout:
                fout.write("%s\n" % inpname)
                fout.write("%s.out\n" % inpname.split(".")[0])
                fout.write("%si\n" % inpname.split(".")[0])
                fout.write("%so\n" % inpname.split(".")[0])
                fout.write("temp\n")
                for element in self.system[0].xyz.specie_labels:
                    fout.write("%s\n" % (element + ".psp8"))
                    #fout.write("%s\n" % (element + ".GGA_PBE-JTH.xml"))
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.chdir("../")

    def system_to_in(self, fout):
        self.system[0].to_in(fout)
        fout.write("\n")
        fout.write("xangst_lastimg\n")
        for atom in self.system[-1].xyz.atoms:
            fout.write("%.9f %.9f %.9f\n" % (atom.x, atom.y, atom.z))
        fout.write("\n")
