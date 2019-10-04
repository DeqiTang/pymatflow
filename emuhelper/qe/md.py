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
        
        self.control.params["calculation"] = 'md'
        self.control.params["outdir"] = "./tmp"
        self.control.params["pseudo_dir"] = "./"
        self.control.params["wf_collect"] = ".true."
        self.system.params["ibrav"] = 0
        self.system.params["nat"] = self.arts.xyz.natom
        self.system.params["ntyp"] = self.arts.xyz.nspecies
        self.system.params["ecutwfc"] = 100
        self.system.params["input_DFT"] = 'PBE'
        self.system.params["occupations"] = 'smearing'
        self.system.params["smearing"] = "gaussian"
        self.system.params["degauss"] = 0.0001
        
        self.ions.params["ion_dynamics"] = 'verlet'
        self.ions.params["ion_temperature"] = 'not_controlled'
        self.ions.params["tempw"] = 300

    def gen_input(self, directory="tmp-md-qe", inpname="molecular-dynamics.in"):
        """
        directory: a place for all the generated files
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        os.system("cp *.UPF %s/" % directory)

        with open(os.path.join(directory, inpname), 'w') as fout:
            self.control.to_in(fout)
            self.system.to_in(fout)
            self.electrons.to_in(fout)
            self.ions.to_in(fout)
            self.arts.to_in(fout)
    
    def run(self, directory="tmp-md-qe", inpname="molecular-dynamics.in", output="molecular-dynamics.out"):
        os.chdir(directory)
        os.system("pw.x < %s | tee %s" % (inpname, output))
        os.chdir("../")

