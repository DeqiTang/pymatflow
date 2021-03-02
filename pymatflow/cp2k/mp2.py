"""
MP2 calculation
"""
import os
import sys
import shutil
import numpy as np


from pymatflow.remote.server import server_handle
from pymatflow.cp2k.cp2k import Cp2k
#from pymatflow.cp2k.base.glob import cp2k_glob
#from pymatflow.cp2k.base.force_eval import cp2k_force_eval
#from pymatflow.cp2k.base.atom import cp2k_atom

"""
"""

class StaticMp2Run(Cp2k):
    """
    Note:
        static_run is the class as an agent for static mp2 type calculation, including
    """
    def __init__(self):
        """
        """
        super().__init__()
        #self.glob = cp2k_glob()
        #self.force_eval = cp2k_force_eval()
        #self.atom = cp2k_atom()

        self.glob.basic_setting(run_type="ENERGY_FORCE")
        self.force_eval.basic_setting()
        self.atom.basic_setting(run_type="MP2")


    def scf_mp2(self, directory="tmp-cp2k-static-mp2", inpname="static-scf-mp2.inp", output="static-scf-mp2.out", runopt="gen", auto=0):
        """
        :param directory:
            directory is and path where the calculation will happen.
        :param inpname:
            input filename for the cp2k
        :param output:
            output filename for the cp2k
        :param force_eval:
            allowing control of FORCE_EVAL/... parameters by user
        :param printout_option:
            a list of integers, controlling the printout of properties, etc.
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

            # using force_eval

            self.force_eval.dft.printout.status = True
            self.force_eval.properties.status = True

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.atom.to_input(fout)

            # gen server job comit file
            self.gen_yh(directory=directory, cmd="$PMF_CP2K", inpname=inpname, output=output)
            # gen pbs server job comit file
            self.gen_pbs(directory=directory, cmd="$PMF_CP2K", inpname=inpname, output=output, jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])

        if runopt == "run" or runopt == "genrun":
           os.chdir(directory)
           os.system("%s $PMF_CP2K -in %s | tee %s" % (self.run_params["mpi"], inpname, output))
           os.chdir("../")
    #
        server_handle(auto=auto, directory=directory, jobfilebase="static-scf-mp2", server=self.run_params["server"])
