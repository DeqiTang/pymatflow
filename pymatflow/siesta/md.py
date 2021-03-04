"""
Molecular Dynamics Calculation
"""
import numpy as np
import sys
import os
import shutil

from pymatflow.remote.server import server_handle
from pymatflow.siesta.siesta import Siesta
#from pymatflow.siesta.base.system import SiestaSystem
#from pymatflow.siesta.base.electrons import SiestaElectrons
#from pymatflow.siesta.base.ions import SiestaIons


class MdRun(Siesta):
    """
    """
    def __init__(self):
        super().__init__()
        #self.system = siesta_system()
        #self.electrons = siesta_electrons()
        #self.ions = siesta_ions()

        self.electrons.basic_setting()
        self.ions.basic_setting(option="md")

    def md(self, directory="tmp-siesta-md", inpname="molecular-dynamics.fdf", output="molecular-dynamics.out",
            runopt="gen", auto=0):
        """
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            shutil.copyfile(self.system.xyz.file, os.path.join(directory, os.path.basename(self.system.xyz.file)))
            for element in self.system.xyz.specie_labels:
                shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))

            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write(self.system.to_string())
                fout.write(self.electrons.to_string())
                fout.write(self.ions.to_string())

            # gen yhbatch script
            self.gen_llhpc(directory=directory, inpname=inpname, output=output, cmd="siesta")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="siesta", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen local bash script
            self.gen_bash(directory=directory, inpname=inpname, output=output, cmd="siesta", mpi=self.run_params["mpi"])


        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            os.system("%s siesta < %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="molecular-dynamics", server=self.run_params["server"])

    #
