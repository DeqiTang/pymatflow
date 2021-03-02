#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import numpy as np
import sys
import shutil

from pymatflow.remote.server import server_handle
from pymatflow.dftbplus.dftbplus import DftbPlus

from pymatflow.vasp.vasp import Vasp

"""
usage:
"""

class MdRun(DftbPlus, Vasp):
    """
    """
    def __init__(self):
        super().__init__()


    def md(self, directory="tmp-dftbplus-aimd", runopt="gen", auto=0):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            #shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.poscar.xyz.file, directory))

            # gen slurm script
            self.gen_llhpc(directory=directory, cmd="$PMF_DFTBPLUS", scriptname="aimd.slurm")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="$PMF_DFTBPLUS", scriptname="aimd.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen local bash script
            self.gen_bash(directory=directory, cmd="$PMF_DFTBPLUS", scriptname="aimd.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="$PMF_DFTBPLUS", scriptname="aimd.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen lsf_sustc script
            self.gen_lsf_sustc(directory=directory, cmd="$PMF_DFTBPLUS", scriptname="aimd.lsf_sustc", jobname=self.run_params["jobname"], np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"], queue=self.run_params["queue"])


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash aimd.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="aimd", server=self.run_params["server"])
    #

