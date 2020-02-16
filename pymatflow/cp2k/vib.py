"""
VIBRATIONAL_ANALYSIS calculation
"""
import numpy as np
import sys
import os
import shutil
import matplotlib.pyplot as plt

from pymatflow.cp2k.cp2k import cp2k

"""
Usage:
"""

class vib_run(cp2k):
    """
    Note:
        vib_run is the calss as an agent for Vibrational Analysis running.
    """
    def __init__(self):
        super().__init__()

        self.glob.basic_setting(run_type="VIBRATIONAL_ANALYSIS")
        self.force_eval.basic_setting()
        # calculation of IR through vib need print dipole moments
        # throught DFT/PRINT/MOMENTS
        #self.force_eval.dft.printout.print_moments()
        self.force_eval.dft.printout.moments.status == True
    

    def vib(self, directory="tmp-cp2k-vib", inpname="vib.inp", output="vib.out", 
            mpi="", runopt="gen",
            jobname="vibrational-analysis", nodes=1, ppn=2):
        """
        :param directory:
            wheere the calculation will happen
        :param inpname:
            input filename for the cp2k
        :param output:
            output filename for the cp2k
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.vibrational_analysis.to_input(fout)
        
            # gen server job comit file
            self.gen_yh(directory=directory, cmd="cp2k.popt", inpname=inpname, output=output)
            # gen pbs server job comit file
            self.gen_pbs(directory=directory, cmd="cp2k.popt", inpname=inpname, output=output, jobname=jobname, nodes=nodes, ppn=ppn)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s cp2k.psmp -in %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")
   # 

