#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import sys
import re
import os
import shutil
import numpy as np
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.qe.base.control import qe_control
from pymatflow.qe.base.system import qe_system
from pymatflow.qe.base.electrons import qe_electrons
from pymatflow.qe.base.arts import qe_arts



"""
Note:
    现在只支持设置ATOMIC_POSITIONS 为crystal类型
    而我喜欢用angstrom, 所以就暂且搁置, 等待以后其支持angstrom
    参考:
    https://atztogo.github.io/phonopy/qe.html
"""




class phonopy_run:
    """
    """
    def __init__(self, xyz_f):
        self.control = qe_control()
        self.system = qe_system()
        self.electrons = qe_electrons()
        self.ions = qe_ions()
        self.arts = qe_arts(xyz_f)

        self.arts.basic_setting(ifstatic=False)
        
        self.system.basic_setting(self.arts)

        self.supercell_n = "1 1 1"
        
    def phonopy(self, directory="tmp-qe-phonopy", pos_inpname="pos.in", head_inpname="head.in",
            mpi="", runopt="gen", control={}, system={}, electrons={}, 
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
            
            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.set_occupations() which uses
            # self.system.set_occupations() to set them, as self.system.set_params() 
            # is suppressed from setting occupations related parameters
            self.set_occupations(system)
            self.control.set_params(control)
            self.system.set_params(system)
            self.electrons.set_params(electrons)
            self.arts.set_kpoints(option=kpoints_option, kpoints_mp=kpoints_mp)

            with open(os.path.join(directory, head_inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.write_kpoints(fout)

            # set up the Phonopy calculation
            os.chdir(directory)
            os.system("cat %s > %s" % (head_inpname, pos_inpname))
            with open(pos_inpname, 'a') as fout:
                self.arts.to_in(fout)
            os.system("phonopy --qe -d --dim='%s' -c %s" % (self.supercell_n, pos_inpname))
            os.system("ls | grep 'supercell-' > pos.data")
            disp_dirs = []
            with open("pos.data", 'r') as fin:
                for line in fin:
                    disp_dirs.append(line.split(".")[0].split("-")[1])
            #with open(head_inpname, 'a') as fout:
            #    self.arts.write_kpoints(fout)
            for disp in disp_dirs:
                os.system("cat %s supercell-%s.in > supercell-%s-full.in" % (head_inpname, disp, disp))
                os.system("rm supercell-%s.in" % disp)
            os.chdir("../")
            # end build the phonopy

            # gen yhbatch script
            #self.gen_yh(directory=directory, inpname=inpname, output=output)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("%s pw.x < %s | tee %s" % (mpi, inpname, output))
            # run the dft
            for disp in disp_dirs:
                os.system("pw.x < supercell-%s-full.in | tee supercell-%s.out" % (disp, disp))

            # analyse the result
            import matplotlib.pyplot as plt

            os.system("phonopy --qe -f supercell-{001..%s}.out" % (disp_dirs[-1]))

            # plot band structure
            with open("band.conf", 'w') as fout:
                fout.write("ATOM_NAME =")
                for element in self.arts.xyz.specie_labels:
                    fout.write(" %s" % element)
                fout.write("\n")
                fout.write("DIM = %s\n" % self.supercell_n)
                fout.write("BAND = 0.5 0.5 0.5 0.0 0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.0\n")
            os.system("phonopy --qe -c %s -p band.conf" % inpname)

            os.chdir("../")
    

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

