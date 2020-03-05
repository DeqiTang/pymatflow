"""
Nudged elastic band calculation
"""
import numpy as np
import sys
import os
import shutil


from pymatflow.remote.server import server_handle
from pymatflow.cp2k.cp2k import cp2k


"""
Note:
    we can check the official neb manual for some information on
    how to run transition state search appropriately.
    usually the inter-image distance between 1~2 Bohr is suggested,
    but I am not sure now whether it is also OK when it is larger
    than 2 Bohr.

    at the beginning, we can use no-CI, and start with a compromised
    scf setting, and restart with a higher precision when it is
    converged(according to the manual, this might increase the
    energy barrier).

"""

class neb_run(cp2k):
    """
    """
    def __init__(self):
        """
        """
        super().__init__()
        #self.glob = cp2k_glob()
        #self.force_eval = cp2k_force_eval()
        #self.motion = cp2k_motion()

        self.glob.basic_setting(run_type="BAND")
        self.force_eval.basic_setting()
        self.motion.set_type("BAND")

    def get_images(self, images):
        """
        :param images:
            ["first.xyz", "intermediate-1.xyz", "intermediate-2.xyz", ..., "last.xyz"]
        """
        self.motion.band.get_images(images)
        self.force_eval.subsys.xyz.get_xyz(images[0])


    def neb(self, directory="tmp-cp2k-neb", inpname="neb.inp", output="neb.out", runopt="gen", auto=0):
        """
        :param directory:
            where the calculation will happen
        :param inpname:
            input filename for the cp2k
        :param output:
            output filename for the cp2k
        :param force_eval:
            allowing control of FORCE_EVAL/... parameters by user
        :param motion:
            allowing control of MOTION/... parameters by user
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            for image in self.motion.band.images:
                shutil.copyfile(image.file, os.path.join(directory, os.path.basename(image.file)))

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.motion.to_input(fout)

            # gen server job comit file
            self.gen_llhpc(directory=directory, inpname=inpname, output=output, cmd="$PMF_CP2K")
            # gen pbs server job comit file
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_CP2K", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_CP2K -in %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="neb", server=self.run_params["server"])
    #
