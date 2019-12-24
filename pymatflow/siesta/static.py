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
        self.properties = siesta_properties(self.system.xyz)
        
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

            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="siesta")

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
           
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="siesta")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s siesta < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")


    def converge_cutoff(self, emin, emax, step, directory="tmp-siesta-cutoff",
            mpi="", runopt="gen", electrons={}, kpoints_mp=[1, 1, 1]):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
        
            for element in self.system.xyz.specie_labels:
                shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))

            self.electrons.kpoints_mp = kpoints_mp
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
            
            # gen yhbatch running script
            with open("converge-cutoff.sub", 'w') as fout:
                fout.write("#!/bin/bash\n")
                for i in range(n_test + 1):
                    meshcutoff = int(emin + i * step)
                    inp_name = "cutoff-%d.fdf" % meshcutoff
                    out_f_name = "cutoff-%d.out" % meshcutoff
                    fout.write("yhrun -N 1 -n 24 siesta < %s > %s\n" % (inp_name, out_f_name))

        if runopt == "run" or runopt == "genrun":
            # run
            os.chdir(directory)
            for i in range(n_test + 1):
                meshcutoff = int(emin + i * step)
                os.system("%s siesta < cutoff-%d.fdf | tee cutoff-%d.out" % (mpi, meshcutoff, meshcutoff))
            os.chdir("../")
    
    def set_spin(self, spin="non-polarized"):
        self.electrons.set_spin(spin)

    def gen_yh(self, inpname, output, directory="tmp-siesta-static", cmd="siesta"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % (cmd, inpname, output))

