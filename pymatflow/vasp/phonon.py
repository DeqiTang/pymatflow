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

class PhononRun(Vasp):
    """
    Reference:
        https://atztogo.github.io/phonopy/vasp-dfpt.html#vasp-dfpt-interface
    """
    def __init__(self):
        super().__init__()

        self.incar.set_runtype(runtype="phonon")

        self.supercell_n = [1, 1, 1]

    def phonon(self, directory="tmp-vasp-phonon", runopt="gen", auto=0):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.poscar.xyz.file, directory))

            ## Construct and run every POSCAR scf
            with open(os.path.join(directory, "POSCAR-original-unitcell"), 'w') as fout:
                self.poscar.to_poscar(fout)

            os.chdir(directory)
            os.system("phonopy -d --dim='%d %d %d' -c POSCAR-original-unitcell" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))
            os.system("cp SPOSCAR POSCAR")
            os.chdir("../")

            #with open(os.path.join(directory, "INCAR"), 'w') as fout:
            #    self.incar.to_incar(fout)

            # gen llhpc script
            self.gen_llhpc(directory=directory, cmd="$PMF_VASP_STD", scriptname="phonon.slurm")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="$PMF_VASP_STD", scriptname="phonon.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen local bash script
            self.gen_bash(directory=directory, cmd="$PMF_VASP_STD", scriptname="phonon.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="$PMF_VASP_STD", scriptname="phonon.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen lsf_sustc script
            self.gen_lsf_sustc(directory=directory, cmd="$PMF_VASP_STD", scriptname="phonon.lsf_sustc", jobname=self.run_params["jobname"], np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cd_cdcloud script
            self.gen_cdcloud(directory=directory, cmd="$PMF_VASP_STD", scriptname="phonon.slurm_cd")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("vasp")
            os.system("bash phonon.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="phonon", server=self.run_params["server"])
    #
