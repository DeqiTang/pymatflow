#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from emuhelper.cp2k.base.glob import cp2k_glob
from emuhelper.cp2k.base.force_eval import cp2k_force_eval

"""
Usage:
"""

class static_run:
    """
    """
    def __init__(self, xyz_f):
        self.glob = cp2k_glob()
        self.force_eval = cp2k_force_eval(xyz_f)
        
        self.glob.basic_setting(run_type="ENERGY_FORCE")
        self.force_eval.basic_setting()


    def scf(self, directory="tmp-cp2k-static", inpname="static-scf.inp", output="static-scf.out", 
            force_eval={}, mpi="", runopt="gen"):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))
            # using force_eval
            self.force_eval.set_params(force_eval)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
    
        if runopt == "run" or runopt == "genrun":
           os.chdir(directory)
           os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
           os.chdir("../")    

    def scf_restart(self, directory="tmp-cp2k-static", inpname="static-scf-restart.inp", output="static-scf-restart.out", 
            force_eval={}, mpi="", runopt="gen"):
        """
        scf_restart continue a scf calculation from previous scf
        or mimic a nscf calculation(there seems no official nscf
        in cp2k) by increasing kpoints from previous scf running
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("scf_restart calculation:\n")
            print("  directory of previous scf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            self.force_eval.dft.scf.params["SCF_GUESS"] = "RESTART"
            # using force_eval
            self.force_eval.set_params(force_eval)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
    
        if runopt == "run" or runopt == "genrun":
           os.chdir(directory)
           os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
           os.chdir("../")    
    
    def converge_cutoff(self, emin, emax, step, rel_cutoff, directory="tmp-cp2k-cutoff", 
            runopt="gen", force_eval={}):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))
        
            n_test = int((emax - emin) / step)
            for i in range(n_test + 1):
                cutoff = int(emin + i * step)
                inpname = "cutoff-%d.inp" % cutoff
                self.force_eval.dft.mgrid.params["CUTOFF"] = cutoff
                self.force_eval.dft.mgrid.params["REL_CUTOFF"] = rel_cutoff
                self.force_eval.set_params(force_eval)
                with open(os.path.join(directory, inpname), 'w') as fout:
                    self.glob.to_input(fout)
                    self.force_eval.to_input(fout)
        # run
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for i in range(n_test + 1):
                cutoff = int(emin + i * step)
                inpname = "cutoff-%d.inp" % cutoff
                output = "cutoff-%d.out" % cutoff
                os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
            os.chdir("../")
            # analyse the result
            os.chdir(directory)
            for i in range(n_test + 1):
                cutoff = int(emin + i * step)
                out_f_name = "cutoff-%d.out" % cutoff
                os.system("cat %s | grep 'Total energy:' >> energy-cutoff.data" % out_f_name)
            cutoffs = [emin + i * step for i in range(n_test + 1)]
            energy = []
            with open("energy-cutoff.data", 'r') as fin:
                for line in fin:
                    energy.append(float(line.split()[2]))
            plt.plot(cutoffs, energy, marker='o')
            plt.title("CUTOFF Converge Test", fontweight='bold', color='red')
            plt.xlabel("CUTOFF (Ry)")
            plt.ylabel("Energy (a.u.)")
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-cutoff.png")
            plt.show()
            os.chdir("../")

    def converge_rel_cutoff(self, emin, emax, step, cutoff, directory="tmp-cp2k-rel-cutoff",
            force_eval={}, runopt="gen"):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))
        
            n_test = int((emax - emin) / step)
            for i in range(n_test + 1):
                rel_cutoff = int(emin + i * step)
                inpname = "rel-cutoff-%d.inp" % rel_cutoff
                self.force_eval.dft.mgrid.params["CUTOFF"] = cutoff
                self.force_eval.dft.mgrid.params["REL_CUTOFF"] = rel_cutoff
                self.force_eval.set_params(force_eval)
                with open(os.path.join(directory, inpname), 'w') as fout:
                    self.glob.to_input(fout)
                    self.force_eval.to_input(fout)
        # run
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for i in range(n_test + 1):
                rel_cutoff = int(emin + i * step)
                inpname = "rel-cutoff-%d.inp" % rel_cutoff
                output = "rel-cutoff-%d.out" % rel_cutoff
                os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
            os.chdir("../")
            # analyse the result
            os.chdir(directory)
            for i in range(n_test + 1):
                rel_cutoff = int(emin + i * step)
                out_f_name = "rel-cutoff-%d.out" % rel_cutoff
                os.system("cat %s | grep 'Total energy:' >> energy-rel-cutoff.data" % out_f_name)
            rel_cutoffs = [emin + i * step for i in range(n_test + 1)]
            energy = []
            with open("energy-rel-cutoff.data", 'r') as fin:
                for line in fin:
                    energy.append(float(line.split()[2]))
            plt.plot(rel_cutoffs, energy, marker='o')
            plt.title("CUTOFF Converge Test", fontweight='bold', color='red')
            plt.xlabel("CUTOFF (Ry)")
            plt.ylabel("Energy (a.u.)")
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-cutoff.png")
            plt.show()
            os.chdir("../")

    def print_electron_density(self):
        self.force_eval.dft.printout.print_electron_density()
    
    def print_bands(self):
        self.force_eval.dft.printout.print_bands()

    def print_pdos(self):
        self.force_eval.dft.printout.print_pdos()
