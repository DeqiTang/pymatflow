#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil
import numpy as np
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.cp2k.base.glob import cp2k_glob
from pymatflow.cp2k.base.force_eval import cp2k_force_eval
#from emuhelper.cp2k.base.atom import cp2k_atom

"""
"""

class static_run:
    """
    Note:
        static_run is the class as an agent for static type calculation, including
        scf, CUTOFF converge test, REL_CUTOFF converge test. and it can control the
        calculation of properties like electronic band structure, projected density
        of states(pdos), electron density, elf, charge difference, etc.
    """
    def __init__(self, xyz_f):
        """
        xyz_f:
            a modified xyz formatted file(the second line specifies the cell of the 
            system).
        TODO: 
            include implement MP2 calculation through CP2K/ATOM
        """
        self.glob = cp2k_glob()
        self.force_eval = cp2k_force_eval(xyz_f)
        #self.atom = cp2k_atom()
        
        self.glob.basic_setting(run_type="ENERGY_FORCE")
        self.force_eval.basic_setting()


    def scf(self, directory="tmp-cp2k-static", inpname="static-scf.inp", output="static-scf.out", 
            force_eval={}, mpi="", runopt="gen", printout_option=[]):
        """
        directory:
            directory is and path where the calculation will happen.
        inpname:
            input filename for the cp2k
        output:
            output filename for the cp2k
        force_eval:
            allowing control of FORCE_EVAL/... parameters by user
        printout_option:
            a list of integers, controlling the printout of properties, etc.
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))
            # using force_eval
            self.force_eval.set_params(force_eval)
            #self.atom.set_params(atom)
            self.force_eval.dft.printout.status = True
            self.force_eval.properties.status = True
            self.printout_option(printout_option)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                #self.atom.to_input(fout)

            # gen server job comit file
            self.gen_yh(cmd="cp2k.popt", inpname=inpname, output=output)
    
        if runopt == "run" or runopt == "genrun":
           os.chdir(directory)
           os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
           os.chdir("../")    

    def scf_restart(self, directory="tmp-cp2k-static", inpname="static-scf-restart.inp", output="static-scf-restart.out", 
            force_eval={}, mpi="", runopt="gen", printout_option=[]):
        """
        Note:
            scf_restart continue a scf calculation from previous scf
            or mimic a nscf calculation(there seems no official nscf
            in cp2k) by increasing kpoints from previous scf running
        directory:
            directory is and path where the calculation will happen.
        inpname:
            input filename for the cp2k
        output:
            output filename for the cp2k
        force_eval:
            allowing control of FORCE_EVAL/... parameters by user
        printout_option:
            a list of integers, controlling the printout of properties, etc.
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
            self.printout_option(printout_option)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
            # gen server job comit file
            #self.gen_yh(cmd="cp2k.popt", inpname=inpname, output=output)   
        if runopt == "run" or runopt == "genrun":
           os.chdir(directory)
           os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
           os.chdir("../")    
    
    def converge_cutoff(self, emin, emax, step, rel_cutoff, directory="tmp-cp2k-cutoff", 
            runopt="gen", force_eval={}):
        """
        Note:
            this function is used to do the converge test for CUTOFF.
            advices on CUTOFF converging:
            we can first check the basis set file, and find the largest
            value among all elements used. and we times it by 4. then
            we can set the converge range convering that value, to find
            the converged CUTOFF
        emin:
            the minimum cutoff of the test range
        emax:
            the maximum cutoff of the test range
        step:
            the step for the converge test
        rel_cutoff:
            for the test of cutoff, rel_cutoff is set to an fixed value
        directory:
            where the converge test happens
        force_eval:
            allowing control of FORCE_EVAL/... parameters by user
        """
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

            # gen yhbatch running script
            with open("converge-cutoff.sub", 'w') as fout:
                fout.write("#!/bin/bash\n")
                for i in range(n_test + 1):
                    cutoff = int(emin + i * step)
                    inpname = "cutoff-%d.inp" % cutoff
                    out_f_name = "cutoff-%d.out" % cutoff
                    fout.write("yhrun -N 1 -n 24 cp2k.popt -in %s | tee %s\n" % (inpname, out_f_name))

        # run
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for i in range(n_test + 1):
                cutoff = int(emin + i * step)
                inpname = "cutoff-%d.inp" % cutoff
                output = "cutoff-%d.out" % cutoff
                os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
            os.chdir("../")

    def converge_rel_cutoff(self, emin, emax, step, cutoff, directory="tmp-cp2k-rel-cutoff",
            force_eval={}, runopt="gen"):
        """
        Note:
            this function is used to do the converge test of REL_CUTOFF.
        emin:
            the minimum rel_cutoff of the test range
        emax:
            the maximum rel_cutoff of the test range
        step:
            the step for the converge test
        cutoff:
            for the test of rel_cutoff, cutoff is set to an fixed value
        directory:
            where the converge test happens
        force_eval:
            allowing control of FORCE_EVAL/... parameters by user
        """
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

            # gen yhbatch running script
            with open("converge-rel-cutoff.sub", 'w') as fout:
                fout.write("#!/bin/bash\n")
                for i in range(n_test + 1):
                    rel_cutoff = int(emin + i * step)
                    inpname = "rel-cutoff-%d.inp" % rel_cutoff
                    out_f_name = "rel-cutoff-%d.out" % rel_cutoff
                    fout.write("yhrun -N 1 -n 24 cp2k.popt -in %s | tee %s\n" % (inpname, out_f_name))

        # run
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for i in range(n_test + 1):
                rel_cutoff = int(emin + i * step)
                inpname = "rel-cutoff-%d.inp" % rel_cutoff
                output = "rel-cutoff-%d.out" % rel_cutoff
                os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
            os.chdir("../")

    def printout_option(self, option=[]):
        """
        Note:
            responsible for the parseing of the printout_option like in self.scf()

        option:
            1: printout pdos
            2: printout band
            3: printout electron densities
            4: printout electron local function(ELF)
            5: printout molecular orbitals
            6: printout molecular orbital cube files
            7: printout mulliken populaltion analysis
            8: printout cubes for generation of STM images
            9: printout cube file with total density(electrons+atomic core)
           10: printout v_hartree_cube
           11: printout v_xc_cube
           12: printout xray_diffraction_spectrum
           13: request a RESP fit of charges.
        """
        if 1 in option:
            self.force_eval.dft.printout.pdos.status = True
        if 2 in option:
            self.force_eval.dft.printout.band_structure.status = True
            self.force_eval.dft.printout.band_structure.set_band(self.force_eval.subsys.xyz)
        if 3 in option:
            self.force_eval.dft.printout.e_density_cube.status = True
        if 4 in option:
            self.force_eval.dft.printout.elf_cube.status = True
        if 5 in option:
            self.force_eval.dft.printout.mo.status = True
        if 6 in option:
            self.force_eval.dft.printout.mo_cubes.status = True
        if 7 in option:
            self.force_eval.dft.printout.mulliken.status = True
        if 8 in option:
            self.force_eval.dft.printout.stm.status = True
        if 9 in option:
            self.force_eval.dft.printout.tot_density_cube.status = True
        if 10 in option:
            self.force_eval.dft.printout.v_hartree_cube.status = True
        if 11 in option:
            self.force_eval.dft.printout.v_xc_cube.status = True
        if 12 in option:
            self.force_eval.dft.printout.xray_diffraction_spectrum.status = True
        if 13 in option:
            self.force_eval.properties.resp.status = True

    def gen_yh(self,inpname, output, directory="tmp-cp2k-static", cmd="cp2k.psmp"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s -in %s | tee %s\n" % (cmd, inpname, output))

