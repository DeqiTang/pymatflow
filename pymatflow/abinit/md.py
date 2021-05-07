"""
Molecular Dynamics calculation
"""
import os
import shutil


from pymatflow.remote.server import server_handle
from pymatflow.abinit.abinit import Abinit

class MdRun(abinit):
    """
    """
    def __init__(self):
        super().__init__()

        self.electrons.basic_setting()

        self.ions.basic_setting(mode="md")

        self.guard.set_queen(queen="md")


    def md(self, directory="tmp-abinit-md", inpname="molecular-dynamics.in", runopt="gen", auto=0):
        #
        self.dataset[0].electrons.set_scf_nscf("scf")

        self.dataset[0].electrons.set_param("tolvrs", 1.0e-10)
        self.dataset[0].electrons.set_param("tolwrf", None)
        self.dataset[0].electrons.set_param("toldff", None)
        self.dataset[0].electrons.set_param("tolrff", None)
        self.dataset[0].electrons.set_param("toldfe", None) #1.0e-6

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
            os.system("cp %s %s/" % (self.dataset[0].system.xyz.file, directory))
            #

            # generate llhpc job submit script
            self.gen_llhpc(directory=directory, script="molecular-dynamics.slurm", cmd="$PMF_ABINIT")
            # generate pbs job submit script
            self.gen_pbs(directory=directory, script="molecular-dynamics.pbs", cmd="$PMF_ABINIT", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # generate local bash job run script
            self.gen_bash(directory=directory, script="molecular-dynamics.sh", cmd="$PMF_ABINIT", mpi=self.run_param["mpi"])
            # generate cdcloud job submit script
            self.gen_cdcloud(directory=directory, script="molecular-dynamics.slurm_cd", cmd="$PMF_ABINIT")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.system("bash %s" % "molecular-dynamics.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="molecular-dynamics", server=self.run_params["server"])

    def analysis(self, directory="tmp-abinit-md", inpname="molecular-dynamics.in"):
        pass
