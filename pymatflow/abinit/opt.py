#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import matplotlib.pyplot as plt

from pymatflow.abinit.abinit import abinit

class opt_run(abinit):
    """
    """
    def __init__(self):
        super().__init__()
        self.input.electrons.basic_setting()
        self.input.ions.basic_setting(mode="opt")

        self.input.guard.set_queen(queen="opt")


    def optimize(self, directory="tmp-abinit-opt", mpi="", runopt="gen",
        jobname="abinit-opt", nodes=1, ppn=32):

        self.input.electrons.set_scf_nscf("scf")

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
            os.system("cp %s %s/" % (self.input.system.xyz.file, directory))

            #
            # generate pbs job submit script
            self.gen_pbs(directory=directory, script="optimization.pbs", cmd="abinit", jobname=jobname, nodes=nodes, ppn=ppn)
            # generate local bash job run script
            self.gen_bash(directory=directory, script="optimization.sh", cmd="abinit", mpi=mpi)


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.system("bash %s" % "optimization.sh")
            os.chdir("../")

    def analysis(self, directory="tmp-abinit-opt", inpname="geometric-optimization.in"):
        pass
