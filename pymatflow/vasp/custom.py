import os
import numpy as np
import sys
import shutil

from pymatflow.remote.server import server_handle
from pymatflow.vasp.vasp import Vasp

"""
usage:
"""

class CustomRun(Vasp):
    """
    Note:
        Custom calculation is not designed for any specific calculation.
        with this type of running, the calculation type is defined by the 
        user provided incar parameters. so this kind of run is flexible.
        for example you can use it to run dfpt calculation.
    """
    def __init__(self):
        super().__init__()

        self.incar.set_runtype(runtype="custom")

    def custom(self, directory="tmp-vasp-custom", runopt="gen", auto=0):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.poscar.xyz.file, directory))

            with open(os.path.join(directory, "INCAR"), 'w') as fout:
                self.incar.to_incar(fout)
            with open(os.path.join(directory, "POSCAR"), 'w') as fout:
                self.poscar.to_poscar(fout)

            # gen slurm script
            self.gen_llhpc(directory=directory, cmd="$PMF_VASP_STD", scriptname="custom.slurm")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="$PMF_VASP_STD", scriptname="custom.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen local bash script
            self.gen_bash(directory=directory, cmd="$PMF_VASP_STD", scriptname="custom.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="$PMF_VASP_STD", scriptname="custom.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen lsf_sustc script
            self.gen_lsf_sustc(directory=directory, cmd="$PMF_VASP_STD", scriptname="custom.lsf_sustc", jobname=self.run_params["jobname"], np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cd_cdcloud script
            self.gen_cdcloud(directory=directory, cmd="$PMF_VASP_STD", scriptname="custom.slurm_cd")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("vasp")
            os.system("bash custom.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="custom", server=self.run_params["server"])
    #
