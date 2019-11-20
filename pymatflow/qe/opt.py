#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.qe.base.control import qe_control
from pymatflow.qe.base.system import qe_system
from pymatflow.qe.base.electrons import qe_electrons
from pymatflow.qe.base.ions import qe_ions
from pymatflow.qe.base.arts import qe_arts


class opt_run:
    """
    structural optimization uses both energies and forces to locate the minima
    along serach directions. usually insufficient scf convergence will lead to
    bad convergence of BFGS algorithm or even to errors. so when doing geometric
    optimization, we better set parameters to get a good scf convergece.
    """
    def __init__(self, xyz_f):
        self.control = qe_control()
        self.system = qe_system()
        self.electrons = qe_electrons()
        self.ions = qe_ions()
        self.arts = qe_arts(xyz_f)

        self.arts.basic_setting(ifstatic=False)
   
        
    def relax(self, directory="tmp-qe-relax", inpname="relax.in", output="relax.out", 
            mpi="", runopt="gen", control={}, system={}, electrons={}, ions={}, 
            kpoints_option="automatic", kpoints_mp=[1, 1, 1, 0, 0, 0]):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.UPF %s/" % directory)
            os.system("cp %s %s/" % (self.arts.xyz.file, directory))
            
            self.set_relax()
            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.set_occupations() which uses
            # self.system.set_occupations() to set them, as self.system.set_params() 
            # is suppressed from setting occupations related parameters
            self.set_occupations(system)
            self.control.set_params(control)
            self.system.set_params(system)
            self.electrons.set_params(electrons)
            self.ions.set_params(ions)
            self.arts.set_kpoints(option=kpoints_option, kpoints_mp=kpoints_mp)

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.ions.to_in(fout)
                self.arts.to_in(fout)
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s pw.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")
    
    def vc_relax(self, directory="tmp-qe-vc-relax", inpname="vc-relax.in", output="vc-relax.out", 
            mpi="", runopt="gen", control={}, system={}, electrons={}, ions={}, 
            kpoints_option="automatic", kpoints_mp=[1, 1, 1, 0, 0, 0]):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.UPF %s/" % directory)
            os.system("cp %s %s/" % (self.arts.xyz.file, directory))

            self.set_vc_relax()
            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.set_occupations() which uses
            # self.system.set_occupations() to set them, as self.system.set_params() 
            # is suppressed from setting occupations related parameters
            self.set_occupations(system)
            self.control.set_params(control)
            self.system.set_params(system)
            self.electrons.set_params(electrons)
            self.ions.set_params(ions)
            self.arts.set_kpoints(option=kpoints_option, kpoints_mp=kpoints_mp)

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.ions.to_in(fout)
                self.arts.to_in(fout)
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output)

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

    def set_occupations(self, system):
        """
            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.system.set_occupations() to 
            # set them, as self.system.set_params() is suppressed from setting
            # occupations related parameters
            # if occupations == None, use default smearing occupation. and 
            # if occupations == "tetrahedra" the value set for smearing and degauss is ignored.
            # if occupations == "smearing", the value of smearing and degauss
            # should be legacy, not None or other illegal values.
        """
        if "occupations" in system:
            if system["occupations"] == None: # user default setting of set_occupations()
                self.system.set_occupations()
            elif system["occupations"] == "tetrahedra":
                self.system.set_occupations(occupations="tetrahedra")
            elif system["occupations"] == "smearing":
                if "smearing" in system and "degauss" in system:
                    self.system.set_occupations(occupations="smearing", smearing=system["smearing"], degauss=system["degauss"])
                elif "smearing" in system:
                    self.system.set_occupations(occupations="smearing", smearing=system["smearing"])
                elif "degauss" in system:
                    self.system.set_occupations(occupations="smearing", degauss=system["degauss"])
                else:
                    self.system.set_occupations(occupations="smearing")
            else:
                pass
    #

    def gen_yh(self, directory, inpname, output):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 pw.x < %s > %s\n" % (inpname, output))

