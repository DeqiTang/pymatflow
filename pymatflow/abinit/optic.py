"""
procedure for Optic calculation:
    1) 
    2) 
    3) 
    4)
    5)
    6)
References:
    * https://docs.abinit.org/tutorial/optic/index.html
    * https://docs.abinit.org/guide/optic/
"""
import os
import sys
import shutil


from pymatflow.remote.server import server_handle
from pymatflow.abinit.abinit import Abinit

from pymatflow.abinit.base.guard import AbinitGuard

from pymatflow.abinit.submit import submit 

class OpticRun(Abinit):
    """
    procedure for Optic calculation:
        1) 
        2) 
        3) 
        4)
        5)
        6)
    Reference:
        * https://docs.abinit.org/tutorial/optic/index.html
        * https://docs.abinit.org/guide/optic/
    Procedure:

    """
    def __init__(self):
        super().__init__()

        self.set_ndtset(6)
    def run(self, directory="matflow-running", runopt="gen", auto=0):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
            os.system("cp %s %s/" % (self.dataset[0].system.xyz.file, directory))

            self.set_pseudos(directory=os.path.abspath(directory))
            """
            self.files.name = "optic-run.files"
            self.files.main_in = "optic-run.in"
            self.files.main_out = "optic-run.out"
            self.files.wavefunc_in = "optic-i"
            self.files.wavefunc_out = "optic-o"
            self.files.tmp = "tmp"
            """

            self.dataset[0].dfpt.status = True
            self.dataset[1].dfpt.status = True
            self.dataset[2].dfpt.status = True
            self.dataset[3].dfpt.status = True
            self.dataset[4].dfpt.status = True
            self.dataset[5].dfpt.status = True
            self.dataset[6].dfpt.status = True

            # overall default dataset 0
            #self.dataset[0].electrons.set_scf_nscf("scf")
            #self.dataset[0].electrons.set_param("nband", 20)
            self.dataset[0].electrons.kpoints.set_param("kptopt", 3)
            self.dataset[0].electrons.set_param("getwfk", 3) # used by 4-6 datasets
            self.dataset[0].electrons.set_param("nbdbuf", 2)

            # 1): SC run with kpoints in the IBZ
            #self.dataset[1].electrons.use_tol(tol="tolvrs", value=1.0e-18)
            self.dataset[1].electrons.use_tol(tol="tolwfr", value=1.0e-22)
            #self.dataset[1].electrons.set_param("nband", 4)
            self.dataset[1].electrons.kpoints.set_param("ngkpt", self.dataset[0].electrons.kpoints.params["ngkpt"].as_val(t=int, dim=1))
            self.dataset[1].misc.set_param("nbdbuf", 0)
            self.dataset[1].misc.set_param("prtden", 1)
            self.dataset[1].electrons.kpoints.set_param("kptopt", 1)
            self.dataset[1].electrons.set_param("getden", 0)
            self.dataset[1].electrons.set_param("getwfk", 0)

            # 2) NSC run with kpoints in the IBZ
            self.dataset[2].electrons.use_tol(tol="tolwfr", value=1.0e-22)
            self.dataset[2].electrons.set_param("iscf", -2)
            self.dataset[2].electrons.set_param("getwfk", 1)
            self.dataset[2].electrons.set_param("getden", 1)
            self.dataset[2].electrons.kpoints.set_param("ngkpt", self.dataset[0].electrons.kpoints.params["ngkpt"].as_val(t=int, dim=1))
            self.dataset[2].electrons.kpoints.set_param("kptopt",1)
            

            # 3) NSC run with large number of bands, and points in the full BZ 
            self.dataset[3].electrons.use_tol(tol="tolwfr", value=1.0e-22)
            self.dataset[3].electrons.kpoints.set_param("ngkpt", self.dataset[0].electrons.kpoints.params["ngkpt"].as_val(t=int, dim=1))
            self.dataset[3].electrons.set_param("iscf", -2)
            self.dataset[3].electrons.set_param("getwfk", 2)
            self.dataset[3].electrons.set_param("getden", 1)
            self.dataset[3].electrons.kpoints.set_param("kptopt",3)

            # 4) ddk response function along axis 1
            #self.dataset[4].electrons.use_tol(tol="tolvrs", value=1.0e-18)
            self.dataset[4].electrons.use_tol(tol="tolwfr", value=1.0e-22)
            self.dataset[4].electrons.kpoints.set_param("ngkpt", self.dataset[0].electrons.kpoints.params["ngkpt"].as_val(t=int, dim=1))
            self.dataset[4].electrons.set_param("iscf", -3)
            self.dataset[4].electrons.set_param("nstep", 1)
            self.dataset[4].misc.set_param("nline", 0)
            self.dataset[4].misc.set_param("prtwf", 3)
            self.dataset[4].dfpt.set_param("nqpt", 1)
            self.dataset[4].dfpt.set_param("rfdir", [1, 0, 0])
            self.dataset[4].dfpt.set_param("rfelfd", 2)
            self.dataset[4].electrons.kpoints.set_param("kptopt",3)
            
            # 5) ddk response function along axis 2
            #self.dataset[5].electrons.use_tol(tol="tolvrs", value=1.0e-18)
            self.dataset[5].electrons.use_tol(tol="tolwfr", value=1.0e-22)
            self.dataset[5].electrons.kpoints.set_param("ngkpt", self.dataset[0].electrons.kpoints.params["ngkpt"].as_val(t=int, dim=1))
            self.dataset[5].electrons.set_param("iscf", -3)
            self.dataset[5].electrons.set_param("nstep", 1)
            self.dataset[5].misc.set_param("nline", 0)
            self.dataset[5].misc.set_param("prtwf", 3)
            self.dataset[5].dfpt.set_param("nqpt", 1)
            self.dataset[5].dfpt.set_param("rfdir", [0, 1, 0])
            self.dataset[5].dfpt.set_param("rfelfd", 2)
            self.dataset[5].electrons.kpoints.set_param("kptopt",3)

            # 6) ddk response function along axis 3
            #self.dataset[6].electrons.use_tol(tol="tolvrs", value=1.0e-18)
            self.dataset[6].electrons.use_tol(tol="tolwfr", value=1.0e-22)
            self.dataset[6].electrons.kpoints.set_param("ngkpt", self.dataset[0].electrons.kpoints.params["ngkpt"].as_val(t=int, dim=1))
            self.dataset[6].electrons.set_param("iscf", -3)
            self.dataset[6].electrons.set_param("nstep", 1)
            self.dataset[6].misc.set_param("nline", 0)
            self.dataset[6].misc.set_param("prtwf", 3)
            self.dataset[6].dfpt.set_param("nqpt", 1)
            self.dataset[6].dfpt.set_param("rfdir", [0, 0, 1])
            self.dataset[6].dfpt.set_param("rfelfd", 2)
            self.dataset[6].electrons.kpoints.set_param("kptopt",3)

        from pymatflow.abinit.optic_namelist import Optic
        submit(abinit=self, directory=directory, prefix="optic-run", runopt=runopt, auto=auto, optic=Optic(system=self.dataset[0].system))