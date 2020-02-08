#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil
import numpy as np
import matplotlib.pyplot as plt

from pymatflow.cp2k.cp2k import cp2k
#from pymatflow.cp2k.base.glob import cp2k_glob
#from pymatflow.cp2k.base.force_eval import cp2k_force_eval
#from pymatflow.cp2k.base.atom import cp2k_atom

"""
"""

class static_mp2_run(cp2k):
    """
    Note:
        static_run is the class as an agent for static mp2 type calculation, including
    """
    def __init__(self):
        """
        xyz_f:
            a modified xyz formatted file(the second line specifies the cell of the 
            system).
        TODO: 
        """
        super().__init__()
        #self.glob = cp2k_glob()
        #self.force_eval = cp2k_force_eval()
        #self.atom = cp2k_atom()
        
        self.glob.basic_setting(run_type="ENERGY_FORCE")
        self.force_eval.basic_setting()
        self.atom.basic_setting(run_type="MP2")


    def scf_mp2(self, directory="tmp-cp2k-static-mp2", inpname="static-scf-mp2.inp", output="static-scf-mp2.out", mpi="", runopt="gen",
            jobname="mp2", nodes=1, ppn=32):
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
           
            self.force_eval.dft.printout.status = True
            self.force_eval.properties.status = True

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.atom.to_input(fout)

            # gen server job comit file
            self.gen_yh(directory=directory, cmd="cp2k.popt", inpname=inpname, output=output)
            # gen pbs server job comit file
            self.gen_pbs(directory=directory, cmd="cp2k.popt", inpname=inpname, output=output, jobname=jobname, nodes=nodes, ppn=ppn)
    
        if runopt == "run" or runopt == "genrun":
           os.chdir(directory)
           os.system("cp2k.psmp -in %s | tee %s" % (inpname, output))
           os.chdir("../")    
    #
