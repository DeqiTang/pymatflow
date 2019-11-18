#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.cp2k.base.glob import cp2k_glob
from pymatflow.cp2k.base.force_eval import cp2k_force_eval
from pymatflow.cp2k.base.motion import cp2k_motion


"""
Usage:
"""

class neb_run:
    """
    """
    def __init__(self, images):
        """
        images:
            ["first.xyz", "intermediate-1.xyz", "intermediate-2.xyz", ..., "last.xyz"]
        """
        self.glob = cp2k_glob()
        self.force_eval = cp2k_force_eval(images[0])
        self.motion = cp2k_motion()
        self.motion.band.get_images(images)

        self.run_type = "BAND"
        self.glob.basic_setting(run_type="BAND")
        self.force_eval.basic_setting()
        self.glob.params["RUN_TYPE"] = "BAND"
        self.motion.set_type("BAND")

    def neb(self, directory="tmp-cp2k-neb", inpname="neb.inp", output="neb.out", 
            mpi="", runopt="gen", force_eval={}, motion={}):
        """
        directory:
            where the calculation will happen
        inpname:
            input filename for the cp2k
        output:
            output filename for the cp2k
        force_eval:
            allowing control of FORCE_EVAL/... parameters by user
        motion:
            allowing control of MOTION/... parameters by user
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            for image in self.motion.band.images:
                shutil.copyfile(image.file, os.path.join(directory, image.file))

            self.force_eval.set_params(force_eval)
            self.motion.set_params(motion)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.motion.to_input(fout)
 
            # gen server job comit file
            self.gen_yh(cmd="cp2k.popt", inpname=inpname, output=output)       

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s cp2k.psmp -in %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")
    
    
    def gen_yh(self,inpname, output, directory="tmp-cp2k-static", cmd="cp2k.psmp"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s -in %s | tee %s\n" % (cmd, inpname, output))
