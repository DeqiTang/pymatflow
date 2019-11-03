#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import matplotlib.pyplot as plt

from emuhelper.abinit.base.electrons import abinit_electrons
from emuhelper.abinit.base.system import abinit_system
from emuhelper.abinit.base.ions import abinit_ions
#from emuhelper.abinit.base.properties import abinit_properties

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

        
    def neb(self, directory="tmp-abinit-neb", inpname="neb.in", mpi="", runopt="gen",
            electrons={}, kpoints={}):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)

            self.electrons.set_scf_nscf("scf")
            self.electrons.set_params(electrons)
            self.electrons.kpoints.set_kpoints(kpoints)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.electrons.to_in(fout)
                self.system_to_in(fout)
                fout.write("nimage 7\n")
                fout.write("imgmov 5\n")
                fout.write("ntimimage 30\n")
                fout.write("tolimg 0.001\n")
                fout.write("dynimage 0 5*1 0\n")
                fout.write("fxcartfactor 1.0\n")

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
