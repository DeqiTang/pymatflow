"""
Linear Response calculation
"""
import numpy as np
import sys
import os
import shutil


from pymatflow.cp2k.cp2k_dev import Cp2k


"""
Usage:
"""

class LrRun(Cp2k):
    """
    Note:
        lr_run is the  class as an agent for Linear Response calculation.
    """
    def __init__(self):
        super().__init__()

        self.set_params({
            "global-run_type": "LINEAR_RESPONSE",
            "force_eval-dft-mgrid-cutoff": 100,
            "force_eval-dft-mgrid-rel_cutoff": 60,
        })


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
            shutil.copyfile(self.xyz.file, os.path.join(directory, os.path.basename(self.xyz.file)))

            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write(self.sections["global"].to_string())
                fout.write(self.sections["force_eval"].to_string())


            # gen server job comit file
            self.gen_yh(cmd="$PMF_CP2K", directory=directory, inpname=inpname, output=output)
            # gen pbs server job comit file
            self.gen_pbs(cmd="$PMF_CP2K", directory=directory, inpname=inpname, output=output, jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])

        if runopt == "run" or runopt == "genrun":
           os.chdir(directory)
           os.system("%s $PMF_CP2K -in %s | tee %s" % (self.run_params["mpi"], inpname, output))
           os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="lr", server=self.run_params["server"])
    #
