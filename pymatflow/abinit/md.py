"""
Molecular Dynamics calculation
"""
import os
import shutil
import matplotlib.pyplot as plt

from pymatflow.remote.server import server_handle
from pymatflow.abinit.abinit import abinit

class md_run(abinit):
    """
    """
    def __init__(self):
        super().__init__()

        self.electrons.basic_setting()

        self.ions.basic_setting(mode="md")

        self.guard.set_queen(queen="md")


    def md(self, directory="tmp-abinit-md", inpname="molecular-dynamics.in", runopt="gen", auto=0):
        #
        self.input.electrons.set_scf_nscf("scf")

        self.input.electrons.params["tolvrs"] = 1.0e-10
        self.input.electrons.params["tolwrf"] = None
        self.input.electrons.params["toldff"] = None
        self.input.electrons.params["tolrff"] = None
        self.input.electrons.params["toldfe"] = None #1.0e-6

        self.files.name = "molecular-dynamics.files"
        self.files.main_in = "molecular-dynamics.in"
        self.files.main_out = "molecular-dynamics.out"
        self.files.wavefunc_in = "md-i"
        self.files.wavefunc_out = "md-o"
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
            self.gen_pbs(directory=directory, script="molecular-dynamics.pbs", cmd="abinit", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # generate local bash job run script
            self.gen_bash(directory=directory, script="molecular-dynamics.sh", cmd="abinit", mpi=self.run_param["mpi"])


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.system("bash %s" % "molecular-dynamics.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="molecular-dynamics", server=self.params["server"])

    def analysis(self, directory="tmp-abinit-md", inpname="molecular-dynamics.in"):
        pass
