"""
calculation of macroscopic polarization through Berry phase approach implemented
in vasp.
from vasp.5.2, we can use LCALPOL and DIPOL to calculate the macroscopi polarization
in a single run, unlike in the past, need three run (IGPAR=1,2,3) to do that.

along with EFIELD_PEAD we might calculate the macroscopic polarization when
finite homogeneous electric field is applied.

if we move the atoms we can also calculate Born effective charge Z* = delta(polarization) / delta(coord)
"""
import os
import sys
import shutil

from pymatflow.remote.server import server_handle
from pymatflow.vasp.vasp import vasp

"""
usage:
"""

class berry_run(vasp):
    """
    """
    def __init__(self):
        super().__init__()

        self.incar.set_runtype(runtype="static")

        self.incar.set_params({
            "LCALCPOL": 'T',
            "DIPOL": [0.25, 0.25, 0.25],
        })

    def berry(self, directory="tmp-vasp-berry", runopt="gen", auto=0):
        """
        directory: a place for all the generated files

        runopt:
            gen    -> generate a new calculation but do not run
            run    -> run a calculation on the previously generated files
            genrun -> generate a calculation and run it

        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.poscar.xyz.file, directory))

            #
            with open(os.path.join(directory, "POSCAR"), 'w') as fout:
                self.poscar.to_poscar(fout)

            # gen yhbatch script
            self.gen_yh(directory=directory, scriptname="polarization-berry.sub", cmd="vasp_std")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="vasp_std", scriptname="polarization-berry.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # gen local bash script
            self.gen_bash(directory=directory, mpi=self.run_params["mpi"], cmd="vasp_std", scriptname="polarization-berry.sh")
            # gen lsf_sz script
            self.gen_lsf_sz(directory=directory, cmd="vasp_std", scriptname="polarization-berry.lsf_sz", np=self.run_params["nodes"]*self.run_params["ppn"], np_per_node=self.run_params["ppn"])

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("vasp")
            os.system("bash polarization-berry.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="polarization-berry", server=self.run_params["server"])

    #
