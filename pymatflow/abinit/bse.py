"""
procedure for BSE calculation:
    1) 
    2) 
    3) 
    4)
    5)
References:
    * https://docs.abinit.org/tutorial/bse/
    * https://docs.abinit.org/theory/bse/
"""
import os
import sys
import shutil


from pymatflow.remote.server import server_handle
from pymatflow.abinit.abinit import Abinit

from pymatflow.abinit.base.guard import AbinitGuard

from pymatflow.abinit.submit import submit 

class BseRun(Abinit):
    """
    procedure for Optic calculation:
        1) 
        2) 
        3) 
        4)
        5)
    Reference:
        * https://docs.abinit.org/tutorial/bse/
        * https://docs.abinit.org/theory/bse/
    Procedure:

    """
    def __init__(self):
        super().__init__()

        self.set_ndtset(5)

    def run(self, directory="matflow-running", runopt="gen", auto=0):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
            os.system("cp %s %s/" % (self.dataset[0].system.xyz.file, directory))

            self.set_pseudos(directory=os.path.abspath(directory))

            #self.dataset[0].dfpt.status = True
            #self.dataset[1].dfpt.status = True
            #self.dataset[2].dfpt.status = True
            #self.dataset[3].dfpt.status = True
            #self.dataset[4].dfpt.status = True
            #self.dataset[5].dfpt.status = True
            #self.dataset[6].dfpt.status = True

            # overall default dataset 0
            #self.dataset[0].electrons.set_scf_nscf("scf")
            #self.dataset[0].electrons.set_param("nband", 20)
            self.dataset[0].electrons.kpoints.set_param("kptopt", 1)
            self.dataset[0].electrons.set_param("nbdbuf", 2)
            self.dataset[0].electrons.kpoints.set_param("nshiftk", 1)
            self.dataset[0].electrons.kpoints.set_param("shiftk", [0.0, 0.0, 0.0]) # Gamma-centered k-mesh
            self.dataset[0].misc.set_param("chkprim", 0)

            # 1): self-consistent calculation
            self.dataset[1].electrons.use_tol(tol="tolvrs", value=1.0e-8)
            self.dataset[1].electrons.kpoints.set_param("ngkpt", self.dataset[0].electrons.kpoints.params["ngkpt"].as_val(t=int, dim=1))
            self.dataset[1].misc.set_param("prtden", 1)

            # 2) definition of parameters for the calculation of the WFK file on the symmetric k-mesh.
            self.dataset[2].electrons.use_tol(tol="tolwfr", value=1.0e-8)
            self.dataset[2].electrons.set_param("iscf", -2) # nscf
            self.dataset[2].electrons.set_param("getden", 1)
            self.dataset[2].electrons.set_param("nbdbuf", 5) # # The last five states are excluded from the converge check to facilitate the convergence
            self.dataset[2].electrons.kpoints.set_param("ngkpt", self.dataset[0].electrons.kpoints.params["ngkpt"].as_val(t=int, dim=1))
            
            # 3) calculation of the WFK file on the shifted k-mesh to break the symmetry.
            self.dataset[3].electrons.use_tol(tol="tolwfr", value=1.0e-8)
            self.dataset[3].electrons.kpoints.set_param("ngkpt", self.dataset[0].electrons.kpoints.params["ngkpt"].as_val(t=int, dim=1))
            self.dataset[3].electrons.set_param("iscf", -2)
            self.dataset[3].electrons.set_param("getden", 1)
            self.dataset[3].electrons.set_param("nbdbuf", 5)
            self.dataset[3].misc.set_param("chksymbreak", 0)
            self.dataset[3].electrons.kpoints.set_param("shiftk", [0.11, 0.21, 0.31])

            # 4) creation of the screening (eps^-1) matrix
            self.dataset[4].electrons.use_tol(tol="tolvrs", value=1.0e-8)
            self.dataset[4].electrons.kpoints.set_param("ngkpt", self.dataset[0].electrons.kpoints.params["ngkpt"].as_val(t=int, dim=1))
            self.dataset[4].electrons.set_param("optdriver", 3)
            self.dataset[4].misc.set_param("gwpara", 2)
            self.dataset[4].misc.set_param("inclvkb", 2)
            self.dataset[4].misc.set_param("awtr", 1)
            self.dataset[4].misc.set_param("symchi", 1)
            self.dataset[4].misc.set_param("getwfk", 2)
            self.dataset[4].misc.set_param("nfreqre", 1) # Only the static limit is needed for standard BSE calculations.
            self.dataset[4].misc.set_param("nfreqim", 0)
            self.dataset[4].electrons.set_param("ecuteps", 6)
            self.dataset[4].electrons.set_param("ecutwfn", 12)

            
            # 5) BSE run with Haydock iterative method (only resonant + W + v)
            #self.dataset[5].electrons.use_tol(tol="tolvrs", value=1.0e-18)
            self.dataset[5].electrons.use_tol(tol="tolwfr", value=1.0e-22)
            self.dataset[5].electrons.kpoints.set_param("ngkpt", self.dataset[0].electrons.kpoints.params["ngkpt"].as_val(t=int, dim=1))
            self.dataset[5].electrons.kpoints.set_param("kptopt", 1)
            self.dataset[5].electrons.kpoints.set_param("nshift", 1)
            self.dataset[5].electrons.kpoints.set_param("shiftk", [0.11, 0.21, 0.31])
            self.dataset[5].misc.set_param("chksymbreak", 0)
            self.dataset[5].electrons.set_param("optdriver", 99) # BS calculation
            self.dataset[5].electrons.set_param("getwfk", 3)
            self.dataset[5].electrons.set_param("getscr", 4)
            self.dataset[5].electrons.set_param("ecutwfn", 12.0)
            self.dataset[5].electrons.set_param("ecuteps", 2.0)
            self.dataset[5].misc.set_param("bs_calctype", 1)
            self.dataset[5].misc.set_param("mbpt_sciss", "0.8 eV")
            self.dataset[5].misc.set_param("bs_exchange_term", 1)
            self.dataset[5].misc.set_param("bs_coulomb_term", 11)
            self.dataset[5].misc.set_param("bs_coupling", 0)
            self.dataset[5].misc.set_param("bs_loband", 2)
            self.dataset[5].misc.set_param("mbpt_sciss", "0 6 0.02 eV")
            self.dataset[5].misc.set_param("bs_algorithm", 2)
            self.dataset[5].misc.set_param("bs_haydock_niter", 200)
            self.dataset[5].misc.set_param("bs_haydock_tol", [0.05, 0])
            self.dataset[5].misc.set_param("zcut", "0.15 eV")
            self.dataset[5].misc.set_param("inclvkb", 2)
            self.dataset[5].misc.set_param("gw_icutcoul", 3)
            self.dataset[5].misc.set_param("istwfk", "*1")

        submit(abinit=self, directory=directory, prefix="bse-run", runopt=runopt, auto=auto)