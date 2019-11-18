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
"""

class opt_run:
    """
    Note:
        opt_run is the calss as an agent for geometric optimization, including GEO_OPT
        and CELL_OPT.
    """
    def __init__(self, xyz_f):
        """
        xyz_f:
            a modified xyz formatted file(the second line specifies the cell of the 
            system).
        TODO: 
        """
        self.glob = cp2k_glob()
        self.force_eval = cp2k_force_eval(xyz_f)
        self.motion = cp2k_motion()
        
        self.run_type = "GEO_OPT" # default is GEO_OPT, can also do CELL_OPT

        self.glob.basic_setting(run_type="ENERGY_FORCE")
        self.force_eval.basic_setting()
        
    def geo_opt(self, directory="tmp-cp2k-geo-opt", inpname="geo-opt.inp", output="geo-opt.out", 
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
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))

            self.set_geo_opt()
            self.force_eval.set_params(force_eval)
            self.motion.set_params(motion)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.motion.to_input(fout)
 
            # gen server job comit file
            self.gen_yh(cmd="cp2k.popt", directory=directory, inpname=inpname, output=output)       

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s cp2k.psmp -in %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")
    
    def cell_opt(self, directory="tmp-cp2k-cell-opt", inpname="cell-opt.inp", output="cell-opt.out", 
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
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))
            
            self.set_cell_opt()
            self.force_eval.set_params(force_eval)
            self.motion.set_params(motion)
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.motion.to_input(fout)
        
            # gen server job comit file
            self.gen_yh(cmd="cp2k.popt", directory=directory, inpname=inpname, output=output)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s cp2k.psmp -in %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def set_geo_opt(self):
        """
        Note:
            set basic parameters for GEO_OPT type running
        """
        self.run_type = "GEO_OPT"
        self.glob.params["RUN_TYPE"] = "GEO_OPT"
        self.motion.set_type("GEO_OPT")
    
    def set_cell_opt(self):
        """
        Note:
            set basic parameters for CELL_OPT type running

        Warning:
            if you are doing CELL_OPT run, you must also enable
            "STRESS_TENSOR" in FORCE_EVAL%STRESS_TENSOR
        """
        self.run_type = "CELL_OPT" 
        self.glob.params["RUN_TYPE"] = "CELL_OPT"
        self.motion.set_type("CELL_OPT")
        if self.force_eval.params["STRESS_TENSOR"] is None:
            self.force_eval.params["STRESS_TENSOR"] = "ANALYTICAL" # NUMERICAL


    def gen_yh(self, inpname, output, directory, cmd="cp2k.psmp"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s -in %s | tee %s\n" % (cmd, inpname, output))
