#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil
import numpy as np

#from pymatflow.cp2k.base.glob import cp2k_glob
#from pymatflow.cp2k.base.force_eval import cp2k_force_eval

from pymatflow.cp2k.cp2k import cp2k

"""
"""

class static_run(cp2k):
    """
    Usage:
        a = static_run()
        a.get_xyz(xxx)
        a.set_params(xxx)
        a.set_printout(xxx)
        a.scf(xxx)
    Note:
        static_run is the class as an agent for static type calculation, including
        scf, CUTOFF converge test, REL_CUTOFF converge test. and it can control the
        calculation of properties like electronic band structure, projected density
        of states(pdos), electron density, elf, charge difference, etc.
    """
    def __init__(self):
        """
        TODO: 
            include implement MP2 calculation through CP2K/ATOM
        """
        super().__init__()
        #self.glob = cp2k_glob()
        #self.force_eval = cp2k_force_eval()
        
        self.glob.basic_setting(run_type="ENERGY_FORCE")
        self.force_eval.basic_setting()

    def get_xyz(self, xyzfile):
        """
        xyz_f:
            a modified xyz formatted file(the second line specifies the cell of the 
            system).
        """
        self.force_eval.subsys.xyz.get_xyz(xyzfile)

    def set_params(self, force_eval={}):
        """
        force_eval:
            allowing control of FORCE_EVAL/... parameters by user
        printout_option:
            a list of integers, controlling the printout of properties, etc.
        """
        # using force_eval
        self.force_eval.set_params(force_eval)            

    def scf(self, directory="tmp-cp2k-static", inpname="static-scf.inp", output="static-scf.out",
            mpi="", runopt="gen"):
        """
        directory:
            directory is and path where the calculation will happen.
        inpname:
            input filename for the cp2k
        output:
            output filename for the cp2k
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)

            # gen server job comit file
            self.gen_yh(cmd="cp2k.popt", inpname=inpname, output=output)
    
        if runopt == "run" or runopt == "genrun":
           os.chdir(directory)
           os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
           os.chdir("../")    
    
    def converge_cutoff(self, emin, emax, step, rel_cutoff, directory="tmp-cp2k-cutoff", runopt="gen"):
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

    def converge_rel_cutoff(self, emin, emax, step, cutoff, directory="tmp-cp2k-rel-cutoff", runopt="gen"):
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


    def converge_kpoints_auto(self, kmin, kmax, step, directory="tmp-cp2k-kpoints-auto", runopt="gen"):
        """
        Note:
            this function is used to do the converge test for KPONTS-AUTO.
            in this test mode, we input kmin, kmax, and step. the code will
            generate continuously the kpoint to test, like when kmin=1, kmax=3
            step=1. the tested kpoints would be: 1x1x1, 2x2x2, 3x3x3.
        kmin:
            the minimum kpoint of the test range
        kmax:
            the maximum kpoint of the test range
        step:
            the step for the converge test
        directory:
            where the converge test happens
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))
        
            n_test = int((kmax - kmin) / step)
            for i in range(n_test + 1):
                kpoint = int(kmin + i * step)
                inpname = "kpoints-%d.inp" % kpoint
                
                self.force_eval.set_params({
                    "DFT-KPOINTS-SCHEME": "MONKHORST-PACK %d %d %d" % (kpoint, kpoint, kpoint)
                    })
                with open(os.path.join(directory, inpname), 'w') as fout:
                    self.glob.to_input(fout)
                    self.force_eval.to_input(fout)

            # gen yhbatch running script
            with open("converge-kpoints.sub", 'w') as fout:
                fout.write("#!/bin/bash\n")
                for i in range(n_test + 1):
                    kpoint = int(kmin + i * step)
                    inpname = "kpoints-%d.inp" % kpoint
                    out_f_name = "kpoints-%d.out" % kpoint
                    fout.write("yhrun -N 1 -n 24 cp2k.popt -in %s | tee %s\n" % (inpname, out_f_name))

        # run
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for i in range(n_test + 1):
                kpoint = int(kmin + i * step)
                inpname = "kpoints-%d.inp" % kpoint
                output = "kpoints-%d.out" % kpoint
                os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
            os.chdir("../")

    def converge_kpoints_manual(self, directory="tmp-cp2k-kpoints-manual", mpi="", runopt="gen",
            kpoints_list=[[1, 1, 1], [2, 2, 2], [3, 3, 3]]):
        """
        Note:
            this function is used to do the converge test for KPOINTS-MANUAL.
            in this mode, we have to specify clearly every k point to test, like
            kpoints_list=[[1, 1, 1], [1, 2, 1], [2, 2, 2]]
        directory:
            where the converge test happens
        force_eval:
            allowing control of FORCE_EVAL/... parameters by user
        kpoints_list:
            kpoints test range
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))
        
            for i in range(len(kpoints_list)):
                inpname = "kpoints-%dx%dx%d.inp" % (kpoints_list[i][0], kpoints_list[i][1], kpoints_list[i][2])
                self.force_eval.set_params({
                    "DFT-KPOINTS-SCHEME": "MONKHORST-PACK %d %d %d" % (kpoints_list[i][0], kpoints_list[i][1], kpoints_list[i][2])
                    })
                with open(os.path.join(directory, inpname), 'w') as fout:
                    self.glob.to_input(fout)
                    self.force_eval.to_input(fout)

            # gen yhbatch running script
            with open("converge-kpoints.sub", 'w') as fout:
                fout.write("#!/bin/bash\n")
                for i in range(len(kpoints_list)):
                    inpname = "kpoints-%dx%dx%d.inp" % (kpoints_list[i][0], kpoints_list[i][1], kpoints_list[i][2])
                    out_f_name = "kpoints-%dx%dx%d.out" % (kpoints_list[i][0], kpoints_list[i][1], kpoints_list[i][2])
                    fout.write("yhrun -N 1 -n 24 cp2k.popt -in %s | tee %s\n" % (inpname, out_f_name))

        # run
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for i in range(len(kpoints_list)):
                inpname = "kpoints-%dx%dx%d.inp" % (kpoints_list[i][0], kpoints_list[i][1], kpoints_list[i][2])
                output = "kpoints-%dx%dx%d.out" % (kpoints_list[i][0], kpoints_list[i][1], kpoints_list[i][2])
                os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
            os.chdir("../")



    def set_printout(self, option=[]):
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
        self.force_eval.dft.printout.status = True
        self.force_eval.properties.status = True

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

    def set_vdw(self, usevdw):
        """
        usevdw: bool
            True or False
        """
        from pymatflow.cp2k.helper import set_vdw
        set_vdw(self.force_eval.dft, usevdw=usevdw)


    def gen_yh(self,inpname, output, directory="tmp-cp2k-static", cmd="cp2k.psmp"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s -in %s | tee %s\n" % (cmd, inpname, output))

