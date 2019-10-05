#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import matplotlib.pyplot as plt

from emuhelper.abinit.base.electrons import abinit_electrons
from emuhelper.abinit.base.system import abinit_system

class static_run:
    """
    """
    def __init__(self, xyz_f):
        self.system = abinit_system(xyz_f)
        self.electrons = abinit_electrons()
        
        self.electrons.params["ecut"] = 50
        self.electrons.params["kptopt"] = 1
        self.electrons.params["ngkpt"] = "1 1 1"
        self.electrons.params["occopt"] = 3  # fermi dirac smearing of occupation
        self.electrons.params["nstep"] = 100
        self.electrons.params["diemac"] = 2.0
        self.electrons.params["toldfe"] = 1.0e-6
        
    def gen_input(self, directory="tmp-abinit-static", inpname="static.in"):
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        os.system("cp *.psp8 %s/" % directory) 

        with open(os.path.join(directory, inpname), 'w') as fout:
            self.electrons.to_in(fout)
            self.system.to_in(fout)

        with open(os.path.join(directory, inpname.split(".")[0]+".files"), 'w') as fout:
            fout.write("%s\n" % inpname)
            fout.write("%s.out\n" % inpname.split(".")[0])
            fout.write("%si\n" % inpname.split(".")[0])
            fout.write("%so\n" % inpname.split(".")[0])
            fout.write("temp\n")
            for element in self.system.xyz.specie_labels:
                fout.write("%s\n" % (element + ".psp8"))
    def run(self, directory="tmp-abinit-static", inpname="static.in"):
        os.chdir(directory)
        os.system("abinit < %s" % inpname.split(".")[0]+".files")
        os.chdir("../")
    
    def analysis(self, directory="tmp-abinit-static", inpname="static.in"):
        pass

    def converge_ecut(self, emin, emax, step, directory="tmp-abinit-ecut"):
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        os.system("cp *.psp8 %s/" % directory)

        os.chdir(directory)

        n_test = int((emax - emin) / step)
        for i in range(n_test + 1):
            cutoff = int(emin + i * step)
            inp_name = "ecut-%d.in" % cutoff
            files_name = "ecut-%d.files" % cutoff
            with open(files_name, 'w') as fout:
                fout.write(inp_name)
                fout.write("\n")
                fout.write("ecut-%d.out\n" % cutoff)
                fout.write("ecut-%di\n" % cutoff)
                fout.write("ecut-%do\n" % cutoff)
                fout.write("temp\n")
                for element in self.system.xyz.specie_labels:
                    fout.write("%s\n" % (element + ".psp8"))
                #
            self.electrons.params["ecut"] = cutoff
            with open(inp_name, 'w') as fout:
                self.electrons.to_in(fout)
                self.system.to_in(fout)
        # run the simulation
        for i in range(n_test + 1):
            cutoff = int(emin + i * step)
            files_name = "ecut-%d.files" % cutoff
            os.system("abinit < %s" % (files_name))

        # analyse the result
        for i in range(n_test + 1):
            cutoff = int(emin + i * step)
            out_f_name = "ecut-%d.out" % cutoff
            os.system("cat %s | grep 'Etotal=' >> energy-ecut.data" % out_f_name)
        ecut = [ emin + i * step for i in range(n_test + 1) ]
        energy = []
        with open("energy-ecut.data", 'r') as fin:
            for line in fin:
                energy.append(float(line.split()[2]))
        #for i in range(len(energy)):
        #    energy[i] = energy[i] - 31
        plt.plot(ecut, energy)
        plt.show()

        os.chdir("../")
        #
