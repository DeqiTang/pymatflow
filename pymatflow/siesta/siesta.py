"""
Overall representation of SIESTA
"""
import numpy as np
import sys
import os
import shutil

from pymatflow.siesta.base.system import SiestaSystem
from pymatflow.siesta.base.electrons import SiestaElectrons
from pymatflow.siesta.base.ions import SiestaIons
from pymatflow.siesta.base.properties import SiestaProperties

class Siesta:
    """
    """
    def __init__(self):
        self.system = SiestaSystem()
        self.electrons = SiestaElectrons()
        self.ions = SiestaIons()
        self.properties = SiestaProperties()

        self.electrons.basic_setting()

        self._initialize()

    def _initialize(self):
        """ initialize the current object, do some default setting
        """
        self.run_params = {}
        self.set_run()

    def get_xyz(self, xyzfile):
        self.system.xyz.get_xyz(xyzfile)
        self.properties.set_xyz(self.system.xyz)

    def set_params(self, params={}, units={}):
        for item in params:
            if item in self.electrons.incharge:
                self.electrons.set_param(item, params[item], unit=units[item] if item in units else None)
            elif item in self.ions.incharge:
                self.ions.set_param(item, params[item], unit=units[item] if item in units else None)
            elif item in self.system.incharge:
                self.system.set_param(item, params[item], unit=units[item] if item in units else None)
            else:
                continue

    def set_kpoints(self, kpoints_mp=[1, 1, 1]):
        self.electrons.kpoints_mp = kpoints_mp

    def set_spin(self, spin="non-polarized"):
        self.electrons.set_spin(spin)

    def set_run(self, mpi="", server="pbs", jobname="siesta", nodes=1, ppn=32, queue=None):
        """ used to set  the parameters controlling the running of the task
        :param mpi: you can specify the mpi command here, it only has effect on native running

        """
        self.run_params["server"] = server
        self.run_params["mpi"] = mpi
        self.run_params["jobname"] = jobname
        self.run_params["nodes"] = nodes
        self.run_params["ppn"] = ppn
        self.run_params["queue"] = queue

    def set_llhpc(self, partition="free", nodes=1, ntask=24, jobname="matflow_job", stdout="slurm.out", stderr="slurm.err"):
        self.run_params["partition"] = partition
        self.run_params["jobname"] = jobname
        self.run_params["nodes"] = nodes
        self.run_params["ntask"] = ntask
        self.run_params["stdout"] = stdout
        self.run_params["stderr"] = stderr

    def gen_llhpc(self, directory, inpname, output, cmd="siesta"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".fdf")[0]+".slurm"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
            fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
            fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
            fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
            fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
            fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
            fout.write("yhrun $PMF_SIESTA\n")


    def gen_yh(self, inpname, output, directory="tmp-siesta-static", cmd="siesta"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".fdf")[0]+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % (cmd, inpname, output))

    def gen_pbs(self, inpname, output, directory, cmd="siesta", jobname="siesta", nodes=1, ppn=32, queue=None):
        """
        generating pbs job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".fdf")[0]+".pbs"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % jobname)
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
            if queue != None:
                fout.write("#PBS -q %s\n" % queue)            
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_SIESTA < %s > %s\n" % (inpname, output))

    def gen_bash(self, inpname, output, directory, cmd="siesta", jobname="siesta", mpi=""):
        """
        generating local bash job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".fdf")[0]+".sh"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % jobname)
            fout.write("%s $PMF_SIESTA < %s > %s\n" % (mpi, inpname, output))
