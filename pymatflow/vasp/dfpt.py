#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil

from pymatflow.remote.server import server_handle
from pymatflow.vasp.vasp import Vasp

"""
usage:
"""

class DfptRun(Vasp):
    """
    Note:

    """
    def __init__(self):
        super().__init__()
        self.incar.set_runtype(runtype="dfpt")


    def dfpt(self, directory="tmp-vasp-dfpt", runopt="gen", mode=0, auto=0):
        """
        directory: a place for all the generated files

        mode:
            0: a brand new calculation, will create a new directory
               and start a brand new dfpt including brand new scf.
            1: there is already a scf calculation, and we will change
               to the existing directory and read the previous data
               and do the dfpt(the original INCAR if exists wil be backup
               as INCAR.old). this will reduce the time needed to
               do the scf calculation.
            dfpt run is based scf ground calculation. but in vasp,
            you can specify scf running and dfpt running in the smame
            INCAR and do them in a single run.
            thus we we will inspect whether there exist a previous scf
            running. if there is, we will not remove the directory but
            just check in and backup the original INCAR to INCAR.old
            and generate the dfpt running INCAR.

        Note:
        """
        if runopt == "gen" or runopt == "genrun":
            if mode == 0:
                if os.path.exists(directory):
                    shutil.rmtree(directory)
                os.mkdir(directory)
                shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
                os.system("cp %s %s/" % (self.poscar.xyz.file, directory))
            elif mode == 1:
                if not os.path.exists(directory):
                    print("================================================\n")
                    print("                 Warning !!!\n")
                    print("================================================\n")
                    print("vasp.dfpt_run.dfpt\n")
                    print("when mode = 1\n")
                    print("there must exists the previous running directory\n")
                    print("but it is not there\n")
                    sys.exit(1)
                #
                if os.path.exists(os.path.join(directory, "INCAR")):
                    os.rename(os.path.join(directory, "INCAR"), os.path.join(directory, "INCAR.old"))
                shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
                os.system("cp %s %s/" % (self.poscar.xyz.file, directory))
            #
            with open(os.path.join(directory, "POSCAR"), 'w') as fout:
                self.poscar.to_poscar(fout)

            # gen llhpc script
            self.gen_llhpc(directory=directory, cmd="$PMF_VASP_STD", scriptname="dfpt.slurm")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="$PMF_VASP_STD", scriptname="dfpt.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen local bash script
            self.gen_bash(directory=directory, cmd="$PMF_VASP_STD", scriptname="dfpt.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="$PMF_VASP_STD", scriptname="dfpt.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen lsf_sustc script
            self.gen_lsf_sustc(directory=directory, cmd="$PMF_VASP_STD", scriptname="dfpt.lsf_sustc", jobname=self.run_params["jobname"], np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cd_cdcloud script
            self.gen_cdcloud(directory=directory, cmd="$PMF_VASP_STD", scriptname="dfpt.slurm_cd")


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("vasp")
            os.system("bash dfpt.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="dfpt", server=self.run_params["server"])


    def set_scf(self, directory="tmp-static-vasp", kpoints=[3, 3, 3]):
        self.kpoints.kpoints = kpoints
        self.incar.set_param("ISTART", 0)
        self.incar.set_param("ICHARG", 2)

        self.incar.to_incar(os.path.join(directory, "INCAR"))
        self.kpoints.to_kpoints(os.path.join(directory, "KPOINTS"))

    #
