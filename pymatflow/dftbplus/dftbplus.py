#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil

from pymatflow.dftbplus.base.geometry import Geometry
from pymatflow.dftbplus.base.hamiltonian import new_hamiltonian
from pymatflow.dftbplus.base.driver import new_driver

from pymatflow.vasp.base.incar import VaspIncar
from pymatflow.vasp.base.poscar import VaspPoscar
from pymatflow.vasp.base.kpoints import VaspKpoints


class DftbPlus:
    """
    """
    def __init__(self):
        self.geometry = Geometry()
        self.hamiltonian = new_hamiltonian()
        self.driver = new_driver()

        self.incar = VaspIncar()
        self.poscar = VaspPoscar()
        self.kpoints = VaspKpoints()

        self._initialize()

    def to_string(self):
        out = ""
        out += self.geometry.to_string()
        out += self.hamiltonian.to_string()
        out += self.driver.to_string()
        return out

    def _initialize(self):
        """ initialize the current object, do some default setting
        """
        self.run_params = {}
        self.set_run()

    def get_xyz(self, xyzfile):
        self.geometry.xyz.get_xyz(xyzfile)

    def set_run(self, mpi="", server="pbs", jobname="cp2k", nodes=1, ppn=32, queue=None):
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

    def gen_llhpc(self, directory, scriptname="dftbplus.sub", cmd="$PMF_DFTBPLUS"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, scriptname), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
            fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
            fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
            fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
            fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
            fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
            fout.write("cat > dftb_in.hsd<<EOF\n")
            fout.write(self.to_string())
            fout.write("EOF\n")
            fout.write("yhrun %s\n" % cmd)


    def gen_yh(self, directory, scriptname="dftbplus.sub", cmd="$PMF_DFTBPLUS"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, scriptname), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("cat > dftb_in.hsd<<EOF\n")
            fout.write(self.to_string())
            fout.write("EOF\n")
            fout.write("yhrun -N 1 -n 24 %s\n" % (cmd))

    def gen_pbs(self, directory, cmd="$PMF_DFTBPLUS", scriptname="dftbplus.pbs", jobname="running", nodes=1, ppn=32, queue=None):
        """
        generating pbs job script for calculation
        """
        with open(os.path.join(directory, scriptname), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % jobname)
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
            if queue != None:
                fout.write("#PBS -q %s\n" % queue)            
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("cat > dftb_in.hsd<<EOF\n")
            fout.write(self.to_string())
            fout.write("EOF\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s \n" % (cmd))            

    def gen_bash(self, directory, mpi="", cmd="$PMF_DFTBPLUS", scriptname="dftbplus.bash"):
        """
        generating bash script for local calculation
        """
        with open(os.path.join(directory, scriptname), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("\n")
            fout.write("cat > dftb_in.hsd<<EOF\n")
            fout.write(self.to_string())
            fout.write("EOF\n")
            fout.write("%s %s\n" % (mpi, cmd))

    def gen_lsf_sz(self, directory, cmd="$PMF_DFTBPLUS", scriptname="dftbplus.lsf_sz", np=24, np_per_node=12, queue="intelY_mid"):
        """
        generating lsf job script for calculation on ShenZhen supercomputer
        """
        with open(os.path.join(directory, scriptname), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("APP_NAME=%s\n" % queue)
            fout.write("NP=%d\n" % np)
            fout.write("NP_PER_NODE=%d\n" % np_per_node)
            fout.write("RUN=\"RAW\"\n")
            fout.write("CURDIR=$PWD\n")
            fout.write("source /home-yg/env/intel-12.1.sh\n")
            fout.write("source /home-yg/env/openmpi-1.6.5-intel.sh\n")
            fout.write("cd $CURDIR\n")
            fout.write("# starting creating ./nodelist\n")
            fout.write("rm -rf $CURDIR/nodelist >& /dev/null\n")
            fout.write("for i in `echo $LSB_HOSTS`\n")
            fout.write("do\n")
            fout.write("  echo \"$i\" >> $CURDIR/nodelist \n")
            fout.write("done\n")
            fout.write("ndoelist=$(cat $CURDIR/nodelist | uniq | awk \'{print $1}\' | tr \'\\n\' \',\')\n")

            fout.write("cat > dftb_in.hsd<<EOF\n")
            fout.write(self.to_string())
            fout.write("EOF\n")
            fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist %s\n" % cmd)

    def gen_lsf_sustc(self, directory, cmd="$PMF_DFTBPLUS", scriptname="dftb_plus.lsf_sustc", jobname="matflow-job", np=24, np_per_node=12, queue="medium"):
        """
        generating lsf job script for calculation on Southern University of Science and Technology supercomputer
        """
        with open(os.path.join(directory, scriptname), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#BSUB -J %s\n" % jobname)
            fout.write("#BSUB -q %s\n" % queue)
            fout.write("#BSUB -n %s\n" % np) #number of total cores
            fout.write("#BSUB -R \"span[ptile=%d]\"\n" % np_per_node)
            fout.write("hostfile=`echo $LSB_DJOB_HOSTFILE`\n")
            fout.write("NP=`cat $hostfile | wc -l`\n")
            fout.write("cd $LS_SUBCWD\n")

            fout.write("cat > dftb_in.hsd<<EOF\n")
            fout.write(self.to_string())
            fout.write("EOF\n")
            fout.write("mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP %s\n" % cmd)


