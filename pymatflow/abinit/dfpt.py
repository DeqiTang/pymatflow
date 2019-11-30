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
from pymatflow.abinit.base.guard import abinit_guard


class dfpt_run():
    """
    procedure for DFPT calculation:
        1) ground-state calculation
        2) DFPT calculation
        3) post-processing
    Reference:
        https://docs.abinit.org/guide/respfn/
        https://docs.abinit.org/tutorial/rf1/
        https://docs.abinit.org/tutorial/rf2/
        https://docs.abinit.org/tutorial/elastic/
        https://docs.abinit.org/tutorial/ffield/
        https://docs.abinit.org/tutorial/nlo/
    Procedure:
        1) first, a self-consistent ground-state computation with the restricted set of k-points in the Irreducible Brillouin Zone (with kptopt=1)

        2) second, the three non-self-consistent response-function computations (one for each direction) of the d/dk perturbation, with the half set of k-points (with kptopt=2, and iscf=-3)

        3) third, all q=0 self-consistent response-function computations of the electric field perturbations and of the atomic displacements, with the half set of k-points (with kptopt=2)

        4) fourth, a non-self-consistent ground-state computation with the set of k+q points, that might be reduced thanks to symmetries (with kptopt=1)

        5) fifth, the self-consistent response-function computations of the atomic displacement perturbations with a q wavevector, with the full set of k-points (with kptopt=3)
    """
    def __init__(self, xyz_f):
        self.system = abinit_system(xyz_f)
        self.electrons = abinit_electrons()
        self.properties = abinit_properties()
        self.dfpt = abinit_dfpt()

        self.electrons.basic_setting()
        #self.dfpt.basic_setting()

        self.guard = abinit_guard(queen="dfpt", electrons=self.electrons, system=self.system, dfpt=self.dfpt)
       
    def run(self, directory="tmp-abinit-static", mpi="", runopt="gen",
            electrons={}, kpoints={}, properties=[]):
        self.nscf_rf_ddk(directory=directory, mpi=mpi, runopt=runopt, electrons=electrons, kpoints=kpoints, properties=properties)
        self.scf_rf_elfd_phon_q0(directory=directory, mpi=mpi, runopt=runopt, electrons=electrons, kpoints=kpoints, properties=properties)
        self.nscf_ground_kq(directory=directory, mpi=mpi, runopt=runopt, electrons=electrons, kpoints=kpoints, properties=properties)
        self.scf_rf_phon_q(directory=directory, mpi=mpi, runopt=runopt, electrons=electrons, kpoints=kpoints, properties=properties)

    def nscf_rf_ddk(self, directory="tmp-abinit-static", inpname="nscf-rf-ddk.in", mpi="", runopt="gen",
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

            #self.electrons.set_scf_nscf("scf")
            self.electrons.set_params(electrons)
            #self.properties.get_option(option=properties)
            self.electrons.params["tolvrs"] = None
            self.electrons.params["toldfe"] = None
            self.electrons.params["toldff"] = None
            self.electrons.params["tolrff"] = None
            self.electrons.params["tolwfr"] = 1.0e-22
            #self.electrons.params["nstep"] = 0
            self.electrons.kpoints.set_params(kpoints)
            self.electrons.params["irdwfk"] = 1
            self.electrons.kpoints.params["kptopt"] = 2
            self.electrons.params["iscf"] = -3
            self.dfpt.params["rfelfd"] = 2
            self.dfpt.params["rfdir"] = [1, 1, 1]
            self.dfpt.params["nqpt"] = 1
            self.dfpt.params["qpt"] = [0, 0, 0]

            #
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.electrons.to_in(fout)
                self.dfpt.to_in(fout)
                #self.properties.to_in(fout)
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

    def scf_rf_elfd_phon_q0(self, directory="tmp-abinit-static", inpname="scf-rf-elfd-phon.in", mpi="", runopt="gen",
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

            #self.electrons.set_scf_nscf("scf")
            self.electrons.set_params(electrons)
            self.electrons.kpoints.set_params(kpoints)
            #self.properties.get_option(option=properties)

            #
            self.electrons.params["nstep"] = None # for scf
            
            self.electrons.params["tolwfr"] = None
            self.electrons.params["tolvrs"] = 1.0e-8
            self.electrons.params["toldfe"] = None
            self.electrons.params["toldff"] = None
            self.electrons.params["tolrff"] = None

            self.electrons.params["irdwfk"] = 1
            self.electrons.params["irdddk"] = 1
            self.electrons.kpoints.params["kptopt"] = 2
            self.electrons.params["iscf"] = None # 7
            self.dfpt.params["rfphon"] = 1
            self.dfpt.params["rfatpol"] = [1, self.system.xyz.natom]
            self.dfpt.params["rfdir"] = [1, 1, 1]
            self.dfpt.params["rfelfd"] = 3
            self.dfpt.params["nqpt"] = 1
            self.dfpt.params["qpt"] = [0, 0, 0]
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.electrons.to_in(fout)
                self.dfpt.to_in(fout)
                #self.properties.to_in(fout)
                self.system.to_in(fout)

            os.system("cp %s/static-scf-output_WFK %s/nscf-rf-ddk-output_WFK" % (directory, directory))
            with open(os.path.join(directory, inpname.split(".")[0]+".files"), 'w') as fout:
                fout.write("%s\n" % inpname)
                fout.write("%s.out\n" % inpname.split(".")[0])
                fout.write("%s-output\n" % "nscf-rf-ddk")
                fout.write("%s-output\n" % inpname.split(".")[0])
                fout.write("temp\n")
                for element in self.system.xyz.specie_labels:
                    fout.write("%s\n" % (element + ".psp8"))
                    #fout.write("%s\n" % (element + ".GGA_PBE-JTH.xml"))
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.chdir("../")


    def nscf_ground_kq(self, directory="tmp-abinit-static", inpname="nscf-ground-kq.in", mpi="", runopt="gen",
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

            #self.electrons.set_scf_nscf("scf")
            #self.electrons.set_params(electrons)
            self.electrons.kpoints.set_params(kpoints)
            #self.properties.get_option(option=properties)
            #
            self.electrons.params["tolwfr"] = 1.0e-22
            self.electrons.params["tolvrs"] = None
            self.electrons.params["toldfe"] = None
            self.electrons.params["toldff"] = None
            self.electrons.params["tolrff"] = None
            self.electrons.params["irdwfk"] = 1
            self.electrons.params["irdddk"] = None #1
            self.electrons.kpoints.params["kptopt"] = 1
            self.electrons.params["iscf"] = -3
            self.dfpt.params["rfphon"] = None #1
            self.dfpt.params["rfatpol"] = None #[1, self.system.xyz.natom]
            self.dfpt.params["rfdir"] = None #[1, 1, 1]
            self.dfpt.params["rfelfd"] = None #3
            self.dfpt.params["nqpt"] = None #1
            self.dfpt.params["qpt"] = None #[0, 0, 0]
            self.dfpt.params["nqpt"] = 1
            self.dfpt.params["qpt"] = [0.5, 0, 0]
            #self.dfpt.params["qptopt"] = 1
            #self.dfpt.params["ngqpt"] = [1, 1, 1]
            #self.guard.check_all()
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

        

    def scf_rf_phon_q(self, directory="tmp-abinit-static", inpname="scf-rf-phon-q.in", mpi="", runopt="gen",
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

            #self.electrons.set_scf_nscf("scf")
            self.electrons.set_params(electrons)
            #self.electrons.kpoints.set_params(kpoints)
            #self.properties.get_option(option=properties)
            #
            self.electrons.params["tolwfr"] = None #1.0e-22
            self.electrons.params["tolvrs"] = 1.0e-8 #None
            self.electrons.params["toldfe"] = None
            self.electrons.params["toldff"] = None
            self.electrons.params["tolrff"] = None
            self.electrons.params["irdwfk"] = 1
            self.electrons.params["irdddk"] = None #1
            self.electrons.params["iscf"] = None # for scf
            self.electrons.kpoints.params["kptopt"] = 3 # for rf calculation at non gamma q point
            self.dfpt.params["rfphon"] = 1
            self.dfpt.params["rfatpol"] = [1, self.system.xyz.natom]
            self.dfpt.params["rfdir"] = [1, 1, 1]
            self.dfpt.params["rfelfd"] = None
            self.dfpt.params["nqpt"] = 1
            self.dfpt.params["qpt"] = [0.5, 0, 0]
            self.dfpt.params["qptopt"] = None  #1
            self.dfpt.params["ngqpt"] = None  #[1, 1, 1] 
            #self.guard.check_all()
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.electrons.to_in(fout)
                self.dfpt.to_in(fout)
                #self.properties.to_in(fout)
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
    #
