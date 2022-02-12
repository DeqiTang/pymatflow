#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil

from pymatflow.remote.server import server_handle
from pymatflow.octopus.octopus import octopus

"""
usage:
"""

class md_run(octopus):
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

            # gen llhpc script
            self.gen_llhpc(directory=directory, cmd="$PMF_VASP_STD", scriptname="md.slurm")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="$PMF_VASP_STD", scriptname="md.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen local bash script
            self.gen_bash(directory=directory, cmd="$PMF_VASP_STD", scriptname="md.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="$PMF_VASP_STD", scriptname="md.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"])

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("vasp")
            os.system("bash md.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="md", server=self.run_params["server"])

    #
