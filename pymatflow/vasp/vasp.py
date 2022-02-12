#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil

from pymatflow.vasp.base.incar import VaspIncar
from pymatflow.vasp.base.poscar import VaspPoscar
from pymatflow.vasp.base.kpoints import VaspKpoints

"""
in the past pymatflow.vasp will generate the INCAR directly,
but that is inconvenient when we want to keep the INCAR
of different kind of calculation. for instance, when you
calculate the band structure you have to continue from
previous scf and nscf calculation, and if you directly
generate the INCAR it will remove the previous INCAR that
is terrible when sometime later you want to check your previous
parameters.
so now pymatflow.vasp will generate the bash script that can
generate the corresponding INCAR for different type of
calculation. and we can check the correspondig bash script
to check the parameter used.
"""

class Vasp:
    """
    """
    def __init__(self):
        self.incar = VaspIncar()
        self.poscar = VaspPoscar()
        self.kpoints = VaspKpoints()

        self._initialize()

    def _initialize(self):
        """ initialize the current object, do some default setting
        """
        self.run_params = {}
        self.set_run()

    def get_xyz(self, xyzfile):
        self.poscar.xyz.get_xyz(xyzfile)

    def set_params(self, params, runtype):
        """
        :param runtype: one onf 'static', 'opt', 'phonopy', 'phonon', 'neb', 'md', 'custom'
        """
        static_allowed = [
            "NWRITE", "PREC", "NCORE", "NPAR", "KPAR", "LPLANE", "NSIM", "ENCUT", "EDIFF", "NELM", "NFREE", "ISMEAR", "SIGMA",
            "IVDW", "LORBIT", "LOPTICS", "CSHIFT", "NEDOS", "ISPIN", "MAGMOM", "LNONCOLLINEAR",
            "LSORBIT","SAXIS", "LMAXMIX", "LHFCALC", "HFSCREEN", "AEXX", "LSUBROT", "ALGO", "IALGO", "ADDGRID",
            "ISYM", "LREAL", "LWAVE", "LCHARG", "LELF", "PSTRESS",
            "AMIX", "BMIX", "LAECHG", "SYMPREC", "NELECT", 
            "KSPACING", "LPARD", "LSEPK", "LSEPB", "NBMOD", "EINT",
            "LVHAR",
        ]
        opt_allowed = [
            "NWRITE", "PREC", "NCORE", "NPAR", "KPAR", "LPLANE", "NSIM", "ENCUT", "EDIFF", "NELM", "NFREE", "ISMEAR", "SIGMA",
            "IVDW", "ISPIN", "MAGMOM", "LNONCOLLINEAR",
            "ALGO", "IALGO", "ADDGRID",
            "ISYM", "LREAL", "LWAVE", "LCHARG", "LAECHG",
            "IBRION", "ISIF", "NSW", "POTIM", "EDIFFG", "PSTRESS",
            "AMIX", "BMIX", "SYMPREC", "NELECT",
            "KSPACING",
            "LVHAR"
        ]
        phonopy_allowed = [
            "NWRITE", "PREC", "NCORE", "NPAR", "KPAR", "LPLANE", "NSIM", "ENCUT", "EDIFF", "NELM", "NFREE", "ISMEAR", "SIGMA",
            "IVDW", "ISPIN", "MAGMOM", "LNONCOLLINEAR",
            "ALGO", "IALGO", "ADDGRID", "LAECHG",
            "ISYM", "LREAL", "LWAVE", "LCHARG",
            "AMIX", "BMIX", "SYMPREC", "NELECT",
            "KSPACING",
        ]
        neb_allowed = [
            "NWRITE", "PREC", "NCORE", "NPAR", "KPAR", "LPLANE", "NSIM", "ENCUT", "EDIFF", "NELM", "NFREE", "ISMEAR", "SIGMA",
            "IVDW", "ISPIN", "MAGMOM", "LNONCOLLINEAR",
            "ALGO", "IALGO", "ADDGRID",
            "ISYM", "LREAL", "LWAVE", "LCHARG",
            "IBRION", "ISIF", "NSW", "POTIM", "EDIFFG",
            "IOPT", "LCLIMB", "SPRING", "IMAGES", "LNEBCELL", "MAXMOVE", "LGLOBAL", 
            "LAUTOSCALE", "INVCURV", "LLINEOPT", "FDSTEP", 
            "AMIX", "BMIX", "SYMPREC", "NELECT",
            "KSPACING",
        ]
        md_allowed = [
            "NWRITE", "PREC", "NCORE", "NPAR", "KPAR", "LPLANE", "NSIM", "ENCUT", "EDIFF", "NELM", "NFREE", "ISMEAR", "SIGMA",
            "IVDW", "ISPIN", "MAGMOM", "LNONCOLLINEAR",
            "ALGO", "IALGO", "ADDGRID",
            "ISYM", "LREAL", "LWAVE", "LCHARG", "LAECHG",
            "IBRION", "ISIF", "NSW", "POTIM", "EDIFFG", "PSTRESS",
            "AMIX", "BMIX", "SYMPREC", "NELECT",
            "KSPACING",
            "MDALGO", "SMASS", "ANDERSON_PROB", "TEBEG", "TEEND",
        ]
        phonon_allowed = [
            "NWRITE", "PREC", "NCORE", "NPAR", "KPAR", "LPLANE", "NSIM", "ENCUT", "EDIFF", "NELM", "NFREE", "ISMEAR", "SIGMA",
            "IVDW", "ISPIN", "MAGMOM", "LNONCOLLINEAR",
            "ALGO", "IALGO", "ADDGRID",
            "ISYM", "LREAL", "LWAVE", "LCHARG", "LAECHG",
            "IBRION", "ISIF", "NSW", "POTIM", "EDIFFG",
            "AMIX", "BMIX", "SYMPREC", "NELECT",
            "KSPACING",
        ]
        if runtype == "static":
            for item in params:
                if item in static_allowed:
                    self.incar.set_params({item: params[item]})
        elif runtype == "opt":
            for item in params:
                if item in opt_allowed:
                    self.incar.set_params({item: params[item]})
        elif runtype == "phonopy":
            for item in params:
                if item in phonopy_allowed:
                    self.incar.set_params({item: params[item]})
        elif runtype == "phonon":
            for item in params:
                if item in phonon_allowed:
                    self.incar.set_params({item: params[item]})
        elif runtype == "neb":
            for item in params:
                if item in neb_allowed:
                    self.incar.set_params({item: params[item]})
        elif runtype == "md":
            for item in params:
                if item in md_allowed:
                    self.incar.set_params({item: params[item]})
        elif runtype == "custom":
            for item in params:
                self.incar.set_params({item: params[item]})
        else:
            pass

    def set_kpoints(self, kpoints_mp=[1, 1, 1, 0, 0, 0], option="automatic",
            kpath=None):
        self.kpoints.set_kpoints(kpoints_mp=kpoints_mp, option=option, kpath=kpath)

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
        #
        self.run_params["partition"] = "default"
        self.run_params["ntask"] = 24
        self.run_params["stdout"] = "xxx.out"
        self.run_params["stderr"] = "xxx.err"

    def set_llhpc(self, partition="free", nodes=1, ntask=24, jobname="matflow_job", stdout="slurm.out", stderr="slurm.err"):
        self.run_params["partition"] = partition
        self.run_params["jobname"] = jobname
        self.run_params["nodes"] = nodes
        self.run_params["ntask"] = ntask
        self.run_params["stdout"] = stdout
        self.run_params["stderr"] = stderr

    def gen_llhpc(self, directory, scriptname="vasp.sub", cmd="vasp_std"):
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
            fout.write("cat > INCAR<<EOF\n")
            self.incar.to_incar(fout)
            fout.write("EOF\n")
            fout.write("cat > KPOINTS<<EOF\n")
            self.kpoints.to_kpoints(fout)
            fout.write("EOF\n")
            fout.write("yhrun %s\n" % cmd)


    def gen_yh(self, directory, scriptname="vasp.sub", cmd="vasp_std"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, scriptname), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("cat > INCAR<<EOF\n")
            self.incar.to_incar(fout)
            fout.write("EOF\n")
            fout.write("cat > KPOINTS<<EOF\n")
            self.kpoints.to_kpoints(fout)
            fout.write("EOF\n")
            fout.write("yhrun -N 1 -n 24 %s\n" % (cmd))

    def gen_pbs(self, directory, cmd="vasp_std", scriptname="vasp.pbs", jobname="vasp", nodes=1, ppn=32, queue=None):
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
            fout.write("cat > INCAR<<EOF\n")
            self.incar.to_incar(fout)
            fout.write("EOF\n")
            fout.write("cat > KPOINTS<<EOF\n")
            self.kpoints.to_kpoints(fout)
            fout.write("EOF\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPI_FABRICS shm:tmi %s \n" % (cmd))
            fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s \n" % (cmd))            

    def gen_bash(self, directory, mpi="", cmd="vasp_std", scriptname="vasp.bash"):
        """
        generating bash script for local calculation
        """
        with open(os.path.join(directory, scriptname), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("\n")
            fout.write("cat > INCAR<<EOF\n")
            self.incar.to_incar(fout)
            fout.write("EOF\n")
            fout.write("cat > KPOINTS<<EOF\n")
            self.kpoints.to_kpoints(fout)
            fout.write("EOF\n")
            fout.write("%s %s\n" % (mpi, cmd))

    def gen_lsf_sz(self, directory, cmd="vasp_std", scriptname="vasp.lsf_sz", np=24, np_per_node=12, queue="intelY_mid"):
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
            fout.write("#VASP=/home-yg/Soft/Vasp5.4/vasp_std\n")
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

            fout.write("cat > INCAR<<EOF\n")
            self.incar.to_incar(fout)
            fout.write("EOF\n")
            fout.write("cat > KPOINTS<<EOF\n")
            self.kpoints.to_kpoints(fout)
            fout.write("EOF\n")
            fout.write("mpirun -np $NP -machinefile $CURDIR/nodelist %s\n" % cmd)

    def gen_lsf_sustc(self, directory, cmd="vasp_std", scriptname="vasp.lsf_sustc", jobname="matflow-job", np=24, np_per_node=12, queue="medium"):
        """
        generating lsf job script for calculation on Southern University of Science and Technology supercomputer
        """
        with open(os.path.join(directory, scriptname), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#BSUB -J %s\n" % jobname)
            fout.write("#BSUB -e %J.err\n")
            fout.write("#BSUB -o %J.out\n")
            fout.write("#BSUB -q %s\n" % queue)
            fout.write("#BSUB -n %s\n" % np) #number of total cores
            fout.write("#BSUB -R \"span[ptile=%d]\"\n" % np_per_node)
            fout.write("hostfile=`echo $LSB_DJOB_HOSTFILE`\n")
            fout.write("NP=`cat $hostfile | wc -l`\n")
            fout.write("cd $LS_SUBCWD\n")

            fout.write("cat > INCAR<<EOF\n")
            self.incar.to_incar(fout)
            fout.write("EOF\n")
            fout.write("cat > KPOINTS<<EOF\n")
            self.kpoints.to_kpoints(fout)
            fout.write("EOF\n")
            fout.write("mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP %s\n" % cmd)


    def set_cdcloud(self, partition="normal", nodes=1, ntask=32, jobname="matflow_job", stdout="slurm.out", stderr="slurm.err"):
        self.run_params["partition"] = partition
        self.run_params["jobname"] = jobname
        self.run_params["nodes"] = nodes
        self.run_params["ntask"] = ntask
        self.run_params["stdout"] = stdout
        self.run_params["stderr"] = stderr

    def gen_cdcloud(self, directory, scriptname="vasp.sub", cmd="vasp_std"):
        """
        generating slurm job script for calculation
        """
        with open(os.path.join(directory, scriptname), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
            fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
            fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
            fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
            fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
            fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
            fout.write("#\n")
            fout.write("export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so\n")
            fout.write("export FORT_BUFFERED=1\n")
            fout.write("cat > INCAR<<EOF\n")
            self.incar.to_incar(fout)
            fout.write("EOF\n")
            fout.write("cat > KPOINTS<<EOF\n")
            self.kpoints.to_kpoints(fout)
            fout.write("EOF\n")
            fout.write("srun --mpi=pmix_v3 %s\n" % cmd)