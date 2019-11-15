#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.siesta.base.system import siesta_system
from pymatflow.siesta.base.electrons import siesta_electrons
from pymatflow.siesta.base.properties import siesta_properties

class static_run:
    """
    """
    def __init__(self, xyz_f):
        self.system = siesta_system(xyz_f)
        self.electrons = siesta_electrons()
        self.properties = siesta_properties()
        
        self.electrons.basic_setting()

                

    def scf(self, directory="tmp-siesta-static", inpname="static-scf.fdf", output="static-scf.out",
            mpi="", runopt="gen", electrons={}, properties=[], kpoints_mp=[1, 1, 1]):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
        
            for element in self.system.xyz.specie_labels:
                shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))
       
            self.electrons.kpoints_mp = kpoints_mp
            self.electrons.set_params(electrons)
            # use self.properties.options to contorl the calculation of properties
            self.properties.options = properties
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.system.to_fdf(fout)
                self.electrons.to_fdf(fout)
                self.properties.to_fdf(fout)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s siesta < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def scf_restart(self, directory="tmp-siesta-static", inpname="static-scf-restart.fdf", output="static-scf-restart.out",
            mpi="", runopt="gen", electrons={}, properties=[], kpoints_mp=[1, 1, 1]):

        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("scf(restart) calculation:\n")
            print("  directory of previous scf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
       
            self.electrons.dm["UseSaveDM"] = "true"
            self.electrons.kpoints_mp = kpoints_mp
            self.electrons.set_params(electrons)
            # use self.properties.option to contorl the calculation of properties
            self.properties.options = properties
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.system.to_fdf(fout)
                self.electrons.to_fdf(fout)
                self.properties.to_fdf(fout)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s siesta < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")


    def analysis(self, directory="tmp-siesta-static", inpname="static.fdf", output="static.out"):
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

    def converge_cutoff(self, emin, emax, step, directory="tmp-siesta-converge-cutoff",
            mpi="", runopt="gen", electrons={}):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
        
            for element in self.system.xyz.specie_labels:
                shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))

            self.electrons.set_params(electrons)
            self.electrons.dm["UseSaveDM"] = "false"

            n_test = int((emax - emin) / step)
            for i in range(n_test + 1):
                meshcutoff = int(emin + i * step)
                self.electrons.params["MeshCutoff"] = meshcutoff
                self.system.label = "siesta-" + str(meshcutoff)
                with open(os.path.join(directory, "cutoff-%d.fdf" % meshcutoff), 'w') as fout:
                    self.system.to_fdf(fout)
                    self.electrons.to_fdf(fout)
                    #self.properties.to_fdf(fout)
        if runopt == "run" or runopt == "genrun":
            # run
            os.chdir(directory)
            for i in range(n_test + 1):
                meshcutoff = int(emin + i * step)
                os.system("%s siesta < cutoff-%d.fdf | tee cutoff-%d.out" % (mpi, meshcutoff, meshcutoff))
            os.chdir("../")
            # analysis
            os.chdir(directory)
            for i in range(n_test + 1):
                meshcutoff = int(emin + i * step)
                out_f_name = "cutoff-%d.out" % meshcutoff
                os.system("cat %s | grep 'Total =' >> energy-cutoff.data" % out_f_name)
            cutoff = [emin + i * step for i in range(n_test + 1) ]
            energy = []
            with open("energy-cutoff.data", 'r') as fin:
                for line in fin:
                    energy.append(float(line.split()[3]))
            plt.plot(cutoff, energy, marker="o")
            plt.title("Energy against MeshCutoff")
            plt.xlabel("MeshCutoff (Ry)")
            plt.ylabel("Energy (Ry)")
            plt.show()
            plt.savefig("energy-cutoff.png")
            os.chdir("../")
    
    def set_spin(self, spin="non-polarized"):
        self.electrons.set_spin(spin)

