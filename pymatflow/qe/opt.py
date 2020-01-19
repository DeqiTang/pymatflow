#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.qe.pwscf import pwscf

class opt_run(pwscf):
    """
    structural optimization uses both energies and forces to locate the minima
    along serach directions. usually insufficient scf convergence will lead to
    bad convergence of BFGS algorithm or even to errors. so when doing geometric
    optimization, we better set parameters to get a good scf convergece.

    when you structure is small, use a large kpoint set, or the optimization
    will not be reliable. if you structure is big enough, a small kpoint set
    will usually suffice the requirement.
    """
    def __init__(self):
        super().__init__()


    def relax(self, directory="tmp-qe-relax", inpname="relax.in", output="relax.out", mpi="", runopt="gen"):
        """
        directory: a place for all the generated files
        """
        self.set_relax()
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            #os.system("cp *.UPF %s/" % directory)
            #os.system("cp %s %s/" % (self.arts.xyz.file, directory))

            # do not copy too many files at the same time or it will be slow
            # so we do not copy all UPF files in the directory but just copy
            # those used in the calculation.
            shutil.copyfile(self.arts.xyz.file, os.path.join(directory, self.arts.xyz.file))
            all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
            for element in self.arts.xyz.specie_labels:
                for upf in all_upfs:
                    if upf.split(".")[0] == element:
                        shutil.copyfile(upf, os.path.join(directory, upf))
                        break
            # 
            
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.ions.to_in(fout)
                self.arts.to_in(fout)
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output)
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s pw.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")
    
    def vc_relax(self, directory="tmp-qe-vc-relax", inpname="vc-relax.in", output="vc-relax.out", mpi="", runopt="gen"):
        """
        directory: a place for all the generated files
        """
        self.set_vc_relax()
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            #os.system("cp *.UPF %s/" % directory)
            #os.system("cp %s %s/" % (self.arts.xyz.file, directory))

            # do not copy too many files at the same time or it will be slow
            # so we do not copy all UPF files in the directory but just copy
            # those used in the calculation.
            shutil.copyfile(self.arts.xyz.file, os.path.join(directory, self.arts.xyz.file))
            all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
            for element in self.arts.xyz.specie_labels:
                for upf in all_upfs:
                    if upf.split(".")[0] == element:
                        shutil.copyfile(upf, os.path.join(directory, upf))
                        break
            # 


            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.ions.to_in(fout)
                self.cell.to_in(fout)
                self.arts.to_in(fout)
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output)
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s pw.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")
        
    def set_relax(self):
        self.control.calculation("relax")
        self.control.basic_setting("relax")
       
        self.system.basic_setting(self.arts)
        self.electrons.basic_setting()
        self.ions.basic_setting()

    def set_vc_relax(self):
        self.control.calculation("vc-relax")
        self.control.basic_setting("vc-relax")
       
        self.system.basic_setting(self.arts)
        self.electrons.basic_setting()
        self.ions.basic_setting()

