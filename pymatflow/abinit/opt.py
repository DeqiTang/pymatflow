#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import matplotlib.pyplot as plt

from pymatflow.remote.server import server_handle
from pymatflow.abinit.abinit import abinit

class opt_run(abinit):
    """
    """
    def __init__(self):
        super().__init__()
        self.dataset[0].electrons.basic_setting()
        self.dataset[0].ions.basic_setting(mode="opt")

        self.dataset[0].guard.set_queen(queen="opt")


    def optimize(self, directory="tmp-abinit-opt", mpi="", runopt="gen", auto=0):

        self.dataset[0].electrons.set_scf_nscf("scf")

        self.files.name = "optimization.files"
        self.files.main_in = "optimization.in"
        self.files.main_out = "optimization.out"
        self.files.wavefunc_in = "optimization-i"
        self.files.wavefunc_out = "optimization-o"
        self.files.tmp = "tmp"

        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
            os.system("cp %s %s/" % (self.dataset[0].system.xyz.file, directory))

            #
            # generate pbs job submit script
            self.gen_pbs(directory=directory, script="optimization.pbs", cmd="abinit", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # generate local bash job run script
            self.gen_bash(directory=directory, script="optimization.sh", cmd="abinit", mpi=self.run_params["mpi"])


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.system("bash %s" % "optimization.sh")
            os.chdir("../")

        server_handle(auto=auto, directory=directory, jobfilebase="optimization", server=self.run_params["server"])

    def analysis(self, directory="tmp-abinit-opt", inpname="geometric-optimization.in"):
        pass
