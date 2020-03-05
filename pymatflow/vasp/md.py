#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil

from pymatflow.remote.server import server_handle
from pymatflow.vasp.vasp import vasp

"""
usage:
"""

class md_run(vasp):
    """
    """
    def __init__(self):
        super().__init__()

        self.incar.set_runtype(runtype="md")

    def md(self, directory="tmp-vasp-md", runopt="gen", ensemble=0, thermostat=0, auto=0):
        """
        directory: a place for all the generated files

        ensemble:
            0: NVE
            1: NVT
            2: NPT
            3: NPH
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.poscar.xyz.file, directory))

            self.incar.set_md(ensemble=ensemble, thermostat=thermostat)

            with open(os.path.join(directory, "POSCAR"), 'w') as fout:
                self.poscar.to_poscar(fout)

            # gen yhbatch script
            self.gen_yh(directory=directory, cmd="vasp", scriptname="md.slurm")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="vasp_std", scriptname="md.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # gen local bash script
            self.gen_bash(directory=directory, cmd="vasp_std", scriptname="md.sh")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("vasp")
            os.system("bash md.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="md", server=self.run_params["server"])

    #
