"""
Linear Response calculation
"""
import numpy as np
import sys
import os
import shutil


from pymatflow.cp2k.cp2k import cp2k
#from pymatflow.cp2k.base.glob import cp2k_glob
#from pymatflow.cp2k.base.force_eval import cp2k_force_eval
#from emuhelper.cp2k.base.atom import cp2k_atom

"""
Usage:
"""

class lr_run(cp2k):
    """
    Note:
        lr_run is the  class as an agent for Linear Response calculation.
    """
    def __init__(self):
        super().__init__()
        #self.glob = cp2k_glob()
        #self.force_eval = cp2k_force_eval()
        #self.atom = cp2k_atom()

        self.glob.basic_setting(run_type="LINEAR_RESPONSE")
        self.force_eval.basic_setting()

    def lr(self, directory="tmp-cp2k-lr", inpname="lr.inp", output="lr.out", runopt="gen", auto=0):
        """
        :param directory:
            where the calculation will happen
        :param inpname:
            inputfile name for the cp2k
        :param output:
            output filename for the cp2k
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                #self.atom.to_input(fout)

            # gen server job comit file
            self.gen_yh(cmd="cp2k.popt", directory=directory, inpname=inpname, output=output)
            # gen pbs server job comit file
            self.gen_pbs(cmd="cp2k.popt", directory=directory, inpname=inpname, output=output, jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])

        if runopt == "run" or runopt == "genrun":
           os.chdir(directory)
           os.system("%s cp2k.psmp -in %s | tee %s" % (self.run_params["mpi"], inpname, output))
           os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="lr", server=self.run_params["server"])
    #
