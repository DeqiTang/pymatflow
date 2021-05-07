"""
Nudged elastic band calculation
"""
import numpy as np
import sys
import os
import shutil


from pymatflow.remote.server import server_handle

from pymatflow.base.xyz import BaseXyz
from pymatflow.cp2k.cp2k_dev import Cp2k


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

class NebRun(Cp2k):
    """
    """
    def __init__(self):
        """
        """
        super().__init__()

        self.set_params({
            "global-run_type": "BAND",
            "force_eval-dft-mgrid-cutoff": 100,
            "force_eval-dft-mgrid-rel_cutoff": 60,
        })
        self.set_section_status({
            "motion-band": True,
        })


        # to generate key output in main output of running
        # if these section are not switched on there are no
        # enough output date of the neb run
        self.set_section_status({
            "motion-band-convergence_info": True,
            "motion-band-program_run_info": True,
            "motion-band-energy": True,
            "motion-band-optimize_band": True,
        })
        self.set_params({
            "motion-band-program_run_info-initial_configuration_info": "TRUE",
        })

        
    def check_neb(self):
        if "OPT_TYPE".lower() not in self.sections["motion"].subsections["band"].subsections["optimize_band"].params or None == self.sections["motion"].subsections["band"].subsections["optimize_band"].params["opt_type"].as_val():
            self.set_section_status({
                "motion-band-optimize_band-diis": True, # use DIIS optimize scheme
                "motion-band-optimize_band-md": False,
            })
        else:
            self.set_section_status({
                "motion-band-optimize_band-diis": False,
                "motion-band-optimize_band-md": True,
            })


    def get_images(self, images):
        """
        :param images:
            ["first.xyz", "intermediate-1.xyz", "intermediate-2.xyz", ..., "last.xyz"]
        """
        for image in images:
            xyz = BaseXyz()
            xyz.get_xyz(image)
            self.images.append(xyz)

        self.xyz.get_xyz(images[0])


    def neb(self, directory="tmp-cp2k-neb", inpname="neb.inp", output="neb.out", runopt="gen", auto=0):
        """
        :param directory:
            where the calculation will happen
        :param inpname:
            input filename for the cp2k
        :param output:
            output filename for the cp2k
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            for image in self.images:
                shutil.copyfile(image.file, os.path.join(directory, os.path.basename(image.file)))

            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write(self.sections["global"].to_string())
                fout.write(self.sections["force_eval"].to_string())
                fout.write(self.sections["motion"].to_string())

            # gen server job comit file
            self.gen_llhpc(directory=directory, inpname=inpname, output=output, cmd="$PMF_CP2K")
            # gen pbs server job comit file
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_CP2K", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cdcloud server job comit file
            self.gen_cdcloud(cmd="$PMF_CP2K", directory=directory, inpname=inpname, output=output)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_CP2K -in %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="neb", server=self.run_params["server"])
    #
