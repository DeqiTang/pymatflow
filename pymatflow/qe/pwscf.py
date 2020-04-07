"""
Overall representation of PWSCF calc
"""
import os
import sys
import shutil

from pymatflow.qe.base.control import qe_control
from pymatflow.qe.base.system import qe_system
from pymatflow.qe.base.electrons import qe_electrons
from pymatflow.qe.base.ions import qe_ions
from pymatflow.qe.base.cell import qe_cell
from pymatflow.qe.base.arts import qe_arts


class pwscf:
    """
    About:
        this class provide a base representation of pwscf
    """
    def __init__(self):
        self.control = qe_control()
        self.system = qe_system()
        self.electrons = qe_electrons()
        self.ions = qe_ions()
        self.cell = qe_cell()
        self.arts = qe_arts()

        #self.control.basic_setting("scf")
        self.electrons.basic_setting()
        self.set_kpoints() # default kpoint setting

        self._initialize()

    def _initialize(self):
        """ initialize the current object, do some default setting
        """
        self.run_params = {}
        self.set_run()

    def get_xyz(self, xyzfile):
        """
        :param xyzfile:
            a modified xyz formatted file(the second line specifies the cell of the
            system).
        """
        self.arts.xyz.get_xyz(xyzfile)
        self.system.basic_setting(self.arts)
        self.arts.basic_setting(ifstatic=True)


    def set_params(self, control={}, system={}, electrons={}, ions={}, cell={}):
        # check if user try to set occupations and smearing and degauss
        # through system. if so, use self.set_occupations() which uses
        # self.system.set_occupations() to set them, as self.system.set_params()
        # is suppressed from setting occupations related parameters
        self.set_occupations(system)
        self.control.set_params(control)
        self.system.set_params(system)
        self.electrons.set_params(electrons)
        self.ions.set_params(ions)
        self.cell.set_params(cell)

    def set_kpoints(self, kpoints_option="automatic", kpoints_mp=[2, 2, 2, 0, 0, 0], crystal_b=None):
        # about the format of crystal_b
        # see corresponding comment for function self.arts.set_kpoints()
        self.arts.set_kpoints(option=kpoints_option, kpoints_mp=kpoints_mp, crystal_b=crystal_b)

    def set_run(self, mpi="", server="pbs", jobname="qe", nodes=1, ppn=32, queue=None):
        """ used to set  the parameters controlling the running of the task
        :param mpi: you can specify the mpi command here, it only has effect on native running

        """
        self.run_params["server"] = server
        self.run_params["mpi"] = mpi
        self.run_params["jobname"] = jobname
        self.run_params["nodes"] = nodes
        self.run_params["ppn"] = ppn
        self.run_params["queue"] = queue


    def set_atomic_forces(self, pressure=None, pressuredir=None):
        self.arts.set_atomic_forces(pressure=pressure, direction=pressuredir)

    def set_occupations(self, system):
        """
            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.system.set_occupations() to
            # set them, as self.system.set_params() is suppressed from setting
            # occupations related parameters
            # if occupations == None, use default smearing occupation. and
            # if occupations == "tetrahedra" the value set for smearing and degauss is ignored.
            # if occupations == "smearing", the value of smearing and degauss
            # should be legacy, not None or other illegal values.
        """
        if "occupations" in system:
            if system["occupations"] == None: # user default setting of set_occupations()
                self.system.set_occupations()
            elif system["occupations"] == "tetrahedra":
                self.system.set_occupations(occupations="tetrahedra")
            elif system["occupations"] == "smearing":
                if "smearing" in system and "degauss" in system:
                    self.system.set_occupations(occupations="smearing", smearing=system["smearing"], degauss=system["degauss"])
                elif "smearing" in system:
                    self.system.set_occupations(occupations="smearing", smearing=system["smearing"])
                elif "degauss" in system:
                    self.system.set_occupations(occupations="smearing", degauss=system["degauss"])
                else:
                    self.system.set_occupations(occupations="smearing")
            elif system["occupations"] == "tetrahedra_lin":
                self.system.set_occupations(occupations="tetrahedra_lin")
            elif system["occupations"] == "tetrahedra_opt":
                self.system.set_occupations(occupations="tetrahedra_opt")
            elif system["occupations"] == "fixed":
                self.system.set_occupations(occupations="fixed")
            elif system["occupations"] == "from_input":
                self.system.set_occupations(occupations="from_input")


    def set_spin(self):
        pass

    def set_llhpc(self, partition="free", nodes=1, ntask=24, jobname="matflow_job", stdout="slurm.out", stderr="slurm.err"):
        self.run_params["partition"] = partition
        self.run_params["jobname"] = jobname
        self.run_params["nodes"] = nodes
        self.run_params["ntask"] = ntask
        self.run_params["stdout"] = stdout
        self.run_params["stderr"] = stderr

    def gen_llhpc(self, inpname, output, directory, cmd="pw.x"):
        """
        generating yhbatch job script for calculation
        better pass in $PMF_PWX
        """
        with open(os.path.join(directory, inpname.split(".in")[0]+".slurm"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
            fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
            fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
            fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
            fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
            fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
            fout.write("yhrun %s < %s > %s\n" % (cmd, inpname, output))


    def gen_yh(self, inpname, output, directory, cmd="pw.x"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".in")[0]+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % (cmd, inpname, output))

    def gen_pbs(self, inpname, output, directory, cmd="pw.x", jobname="pwscf", nodes=1, ppn=32, queue=None):
        """
        generating pbs job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".in")[0]+".pbs"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % jobname)
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
            if queue != None:
                fout.write("#PBS -q %s\n" % queue)            
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))
