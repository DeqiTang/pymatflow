#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from emuhelper.qe.base.control import qe_control
from emuhelper.qe.base.system import qe_system
from emuhelper.qe.base.electrons import qe_electrons
from emuhelper.qe.base.ions import qe_ions
from emuhelper.qe.base.arts import qe_arts


class md_run:
    """
    """
    def __init__(self, xyz_f):
        self.control = qe_control()
        self.system = qe_system()
        self.electrons = qe_electrons()
        self.ions = qe_ions()
        self.arts = qe_arts(xyz_f)
        
        
    def md(self, directory="tmp-qe-md", inpname="md.in", output="md.out", mpi="", runopt="gen",
            control={}, system={}, electrons={}, ions={}, kpoints_mp=[1, 1, 1, 0, 0, 0]):
        """
        directory: a place for all the generated files
        """
        if runopt ==  "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.UPF %s/" % directory)
            
            self.set_md()
            self.control.set_params(control)
            self.system.set_params(system)
            self.electrons.set_params(electrons)
            self.ions.set_params(ions)
            self.arts.set_kpoints(kpoints_mp)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.ions.to_in(fout)
                self.arts.to_in(fout)
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s pw.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def vc_md(self, directory="tmp-qe-vc-md", inpname="vc-md.in", output="vc-md.out", mpi="", runopt="gen", 
            control={}, system={}, electrons={}, ions={}, kpoints_mp=[1, 1, 1, 0, 0, 0]):
        """
        directory: a place for all the generated files
        """
        if runopt ==  "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.UPF %s/" % directory)
            
            self.set_vc_md()
            self.control.set_params(control)
            self.system.set_params(system)
            self.electrons.set_params(electrons)
            self.ions.set_params(ions)
            self.arts.set_kpoints(kpoints_mp)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.ions.to_in(fout)
                self.arts.to_in(fout)
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s pw.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def set_md(self):
        self.control.calculation('md')
        self.control.basic_setting("md")
        
        self.system.basic_setting(self.arts)
        self.electrons.basic_setting()
        self.ions.basic_setting('md') 

    def set_vc_md(self):
        self.control.calculation('vc-md')
        self.control.basic_setting("vc-md")
        
        self.system.basic_setting(self.arts)
        self.electrons.basic_setting()
        self.ions.basic_setting('vc-md') 

