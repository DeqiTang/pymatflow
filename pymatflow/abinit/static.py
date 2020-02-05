#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import matplotlib.pyplot as plt

from pymatflow.abinit.abinit import abinit
#from pymatflow.abinit.base.electrons import abinit_electrons
#from pymatflow.abinit.base.system import abinit_system
#from pymatflow.abinit.base.properties import abinit_properties

class static_run(abinit):
    """
    GOAL: support for both single dataset and multi-dataset mode in abinit,
          currently, only for single dataset mode
    """
    def __init__(self):
        super().__init__()
        #self.system = abinit_system()
        #self.electrons = abinit_electrons()
        #self.properties = abinit_properties()
        self.guard.set_queen(queen="static", electrons=self.electrons, system=self.system)

        self.electrons.basic_setting()

        
    def scf(self, directory="tmp-abinit-static", inpname="static-scf.in", mpi="", runopt="gen"):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
            os.system("cp %s %s/" % (self.system.xyz.file, directory))

            self.electrons.set_scf_nscf("scf")
            #
            self.guard.check_all()
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.electrons.to_in(fout)
                self.properties.to_in(fout)
                self.system.to_in(fout)

            with open(os.path.join(directory, inpname.split(".")[0]+".files"), 'w') as fout:
                fout.write("%s\n" % inpname)
                fout.write("%s.out\n" % inpname.split(".")[0])
                fout.write("%s-input\n" % inpname.split(".")[0])
                fout.write("%s-output\n" % inpname.split(".")[0])
                fout.write("temp\n")
                for element in self.system.xyz.specie_labels:
                    fout.write("%s\n" % (element + ".psp8"))
                    #fout.write("%s\n" % (element + ".GGA_PBE-JTH.xml"))
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.chdir("../")
 
    def nscf(self, directory="tmp-abinit-static", inpname="static-nscf.in", mpi="", runopt="gen"):
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("non-scf calculation:\n")
            print("  directory of previous scf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":

            self.electrons.set_scf_nscf("nscf")
            self.electrons.params["irdwfk"] = 1
            self.electrons.params["irdden"] = 1
            # dos
            self.electrons.params["nband"] = 10
            self.electrons.params["dosdeltae"] = 0.00005
            self.electrons.params["occopt"] = 7
            self.electrons.params["tsmear"] = 0.0001
            # end dos
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.electrons.to_in(fout)
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
 

    def band(self, directory="tmp-abinit-static", inpname="static-band.in", mpi="", runopt="gen"):
        """
        we can use abiopen.py static-band-output_GSR.nc --expose -sns=talk to view the band structure.
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("band structure calculation:\n")
            print("  directory of previous scf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":

            self.electrons.params["iscf"] = -2
            self.electrons.params["nband"] = 8
            self.electrons.params["tolwfr"] = 1.0e-12 # when kptopt < 0 namely band structure calculatin, we can only use tolwfr
            self.electrons.params["tolvrs"] = None
            self.electrons.params["toldfe"] = None
            #self.electrons.params["irdden"] = 1 # actually irdden will be 1 by default if iscf < 0
            
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.electrons.to_in(fout)
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

    def converge_ecut(self, emin, emax, step, directory="tmp-abinit-ecut", mpi="", runopt="gen"):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp %s %s/" % (self.system.xyz.file, directory))
   
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
                    fout.write("ecut-%d-input\n" % cutoff)
                    fout.write("ecut-%d-output\n" % cutoff)
                    fout.write("temp\n")
                    for element in self.system.xyz.specie_labels:
                        fout.write("%s\n" % (element + ".psp8"))
                #
                self.electrons.params["ecut"] = cutoff
                with open(inp_name, 'w') as fout:
                    self.electrons.to_in(fout)
                    self.system.to_in(fout)
            os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            for i in range(n_test + 1):
                cutoff = int(emin + i * step)
                files_name = "ecut-%d.files" % cutoff
                os.system("abinit < %s" % (files_name))
            os.chdir("../")

            # analyse the result
            os.chdir(directory)
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
            plt.title("Ecut converge test")
            plt.xlabel("Ecut (Hartree)")
            plt.ylabel("Energy")
            plt.savefig("energy-ecut.png")
            plt.show()

            os.chdir("../")
        #
    #
