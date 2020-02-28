"""
# ==============================================================================
# pymatflow.abinit.dfpt:
# the general control for dfpt type running
# ==============================================================================
procedure for DFPT calculation:
    1) ground-state calculation
    2) DFPT calculation
    3) post-processing
Reference:
    https://docs.abinit.org/guide/respfn/
    https://docs.abinit.org/tutorial/rf1/
    https://docs.abinit.org/tutorial/rf2/
    https://docs.abinit.org/tutorial/elastic/
    https://docs.abinit.org/tutorial/ffield/
    https://docs.abinit.org/tutorial/nlo/
Procedure:
    1) first, a self-consistent ground-state computation with the restricted set of k-points in the Irreducible Brillouin Zone (with kptopt=1)

    2) second, the three non-self-consistent response-function computations (one for each direction) of the d/dk perturbation, with the half set of k-points (with kptopt=2, and iscf=-3)

    3) third, all q=0 self-consistent response-function computations of the electric field perturbations and of the atomic displacements, with the half set of k-points (with kptopt=2)

    4) fourth, a non-self-consistent ground-state computation with the set of k+q points, that might be reduced thanks to symmetries (with kptopt=1)

    5) fifth, the self-consistent response-function computations of the atomic displacement perturbations with a q wavevector, with the full set of k-points (with kptopt=3)

"""
import os
import sys
import shutil
import matplotlib.pyplot as plt

from pymatflow.remote.server import server_handle
from pymatflow.abinit.abinit import abinit

from pymatflow.abinit.base.guard import abinit_guard


class dfpt_elastic_piezo(abinit):
    """
    procedure for DFPT calculation:
        1) ground-state calculation
        2) DFPT calculation
        3) post-processing
    Reference:
        https://docs.abinit.org/guide/respfn/
        https://docs.abinit.org/tutorial/rf1/
        https://docs.abinit.org/tutorial/rf2/
        https://docs.abinit.org/tutorial/elastic/
        https://docs.abinit.org/tutorial/ffield/
        https://docs.abinit.org/tutorial/nlo/
    Procedure:

    """
    def __init__(self):
        super().__init__()
        #self.dataset[0].guard.set_queen(queen="dfpt")

        self.set_ndtset(3)

    def run(self, directory="tmp-abinit-dfpt-elastic-piezo", runopt="gen", auto=0):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
            os.system("cp %s %s/" % (self.dataset[0].system.xyz.file, directory))

            self.files.name = "dfpt-elastic-piezo.files"
            self.files.main_in = "dfpt-elastic-piezo.in"
            self.files.main_out = "dfpt-elastic-piezo.out"
            self.files.wavefunc_in = "dfpt-elastic-piezo-i"
            self.files.wavefunc_out = "dfpt-elastic-piezo-o"
            self.files.tmp = "tmp"

            self.dataset[0].dfpt.status = True
            self.dataset[1].dfpt.status = True
            self.dataset[2].dfpt.status = True
            self.dataset[3].dfpt.status = True

            # overall defaultg dataset 0
            #self.dataset[0].electrons.set_scf_nscf("scf")
            # 1) ground state scf calculation
            self.dataset[1].electrons.use_tol(tol="tolvrs", value=1.0e-10)
            self.dataset[1].electrons.kpoints.params["ngkpt"] = self.dataset[0].electrons.kpoints.params["ngkpt"]

            # 2) calculate ddk wave functions (iscf = -3, rfelfd =2, qpt = 0 0 0, rfdir = 1 1 1)
            # note here rfelfd = 2 means we only calculate perturbation to dk and
            # perturbation to electric field is not accounted! see https://docs.abinit.org/tutorial/elastic/
            self.dataset[2].electrons.use_tol(tol="tolwfr", value=1.0e-22)
            self.dataset[2].electrons.params["getwfk"] = -1
            self.dataset[2].electrons.kpoints.params["kptopt"] = 2
            self.dataset[2].electrons.kpoints.params["ngkpt"] = self.dataset[0].electrons.kpoints.params["ngkpt"]
            self.dataset[2].electrons.params["iscf"] = -3
            self.dataset[2].dfpt.params["rfelfd"] = 2
            self.dataset[2].dfpt.params["rfdir"] = [1, 1, 1]
            self.dataset[2].dfpt.params["nqpt"] = 1
            self.dataset[2].dfpt.params["qpt"] = [0, 0, 0]


            # 3) calculate 2DTE for elastic and piezoelectric tensors (rfphon = 1, rfatpol, rfdir = 1 1 1, rfstrs = 3)
            self.dataset[3].electrons.use_tol(tol="tolvrs", value=1.0e-10)
            self.dataset[3].electrons.kpoints.params["kptopt"] = 2
            self.dataset[3].electrons.kpoints.params["ngkpt"] = self.dataset[0].electrons.kpoints.params["ngkpt"]
            self.dataset[3].electrons.params["getddk"] = -1
            self.dataset[3].electrons.params["getwfk"] = -2
            self.dataset[3].dfpt.params["nqpt"] = 1
            self.dataset[3].dfpt.params["qpt"] = [0, 0, 0]
            self.dataset[3].dfpt.params["rfphon"] = 1
            self.dataset[3].dfpt.params["rfelfd"] = None
            self.dataset[3].dfpt.params["rfatpol"] = [1, self.dataset[0].system.xyz.natom]
            self.dataset[3].dfpt.params["rfdir"] = [1, 1, 1]
            self.dataset[3].dfpt.params["rfstrs"] = 3


            # pbs jobsubmit script
            with open(os.path.join(directory, "dfpt-elastic-piezo.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                fout.write("cat > %s<<EOF\n" % self.files.main_in)
                fout.write(self.to_string())
                fout.write("EOF\n")
                fout.write("cat > %s<<EOF\n" % self.files.name)
                fout.write(self.files.to_string(system=self.dataset[0].system))
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s\n" % ("abinit", self.files.name))

                # use anaddb to analyse DDB file generated previously
                #with open(os.path.join(directory, "anaddb.in"), "w") as fout:
                fout.write("cat > %s <<EOF\n" % "anaddb.in")
                fout.write("elaflag 3\n")
                fout.write("piezoflag 3\n")
                fout.write("instrflag 1\n")
                fout.write("chneut 1\n")
                fout.write("EOF\n")
                #with open(os.path.join(directory, "anaddb.files"), 'w') as fout:
                fout.write("cat > %s <<EOF\n" % "anaddb.files")
                fout.write("anaddb.in\n")
                fout.write("anaddb.out\n")
                fout.write("%s_DS3_DDB\n" % self.files.wavefunc_out)
                fout.write("dummy_moldyn\n")
                fout.write("dummy_GKK\n")
                fout.write("dummy_epout\n")
                fout.write("dummy_ddk\n")
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s\n" % ("anaddb", "anaddb.files"))

            # pbs local bash script
            with open(os.path.join(directory, "dfpt-elastic-piezo.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n")

                fout.write("cat > %s<<EOF\n" % self.files.main_in)
                fout.write(self.to_string())
                fout.write("EOF\n")
                fout.write("cat > %s<<EOF\n" % self.files.name)
                fout.write(self.files.to_string(system=self.dataset[0].system))
                fout.write("EOF\n")
                fout.write("%s %s < %s\n" % (self.run_params["mpi"], "abinit", self.files.name))

                # use anaddb to analyse DDB file generated previously
                #with open(os.path.join(directory, "anaddb.in"), "w") as fout:
                fout.write("cat > %s <<EOF\n" % "anaddb.in")
                fout.write("elaflag 3\n")
                fout.write("piezoflag 3\n")
                fout.write("instrflag 1\n")
                fout.write("chneut 1\n")
                fout.write("EOF\n")
                #with open(os.path.join(directory, "anaddb.files"), 'w') as fout:
                fout.write("cat > %s <<EOF\n" % "anaddb.files")
                fout.write("anaddb.in\n")
                fout.write("anaddb.out\n")
                fout.write("%s_DS3_DDB\n" % self.files.wavefunc_out)
                fout.write("dummy_moldyn\n")
                fout.write("dummy_GKK\n")
                fout.write("dummy_epout\n")
                fout.write("dummy_ddk\n")
                fout.write("EOF\n")
                fout.write("%s %s < %s\n" % (self.run_params["mpi"], "anaddb", "anaddb.files"))

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash dfpt-elastic-piezo.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="dfpt-elastic-piezo", server=self.run_params["server"])




class dfpt_phonon(abinit):
    """
    procedure for DFPT calculation:
        1) ground-state calculation
        2) DFPT calculation
        3) post-processing
    Reference:
        https://docs.abinit.org/guide/respfn/
        https://docs.abinit.org/tutorial/rf1/
        https://docs.abinit.org/tutorial/rf2/
        https://docs.abinit.org/tutorial/elastic/
        https://docs.abinit.org/tutorial/ffield/
        https://docs.abinit.org/tutorial/nlo/
    Procedure:
    """
    def __init__(self):
        super().__init__()
        #self.dataset[0].electrons.basic_setting()

    def get_qpath(self, qpath):
        self.qpath = qpath

        # decide howmany dataset to use based on qpath
        nq = 0
        specialk = []
        for kpoint in self.qpath:
            if kpoint[3] not in specialk:
                nq += 1
                specialk.append(kpoint[3])
            else:
                pass
        if 'GAMMA' in specialk:
            # then we can calculate for q = GAMMA during esponse function
            # calculation of Q=0 phonons and electric field perturbation
            # and no need to use a dataset to calculated phonon on q = GAMMA
            self.set_ndtset(nq+2)
        else:
            self.set_ndtset(nq+3)

    def run(self, directory="tmp-abinit-dfpt-phonon", runopt="gen", auto=0):

        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
            os.system("cp %s %s/" % (self.dataset[0].system.xyz.file, directory))

            self.files.name = "dfpt-phonon.files"
            self.files.main_in = "dfpt-phonon.in"
            self.files.main_out = "dfpt-phonon.out"
            self.files.wavefunc_in = "dfpt-phonon-i"
            self.files.wavefunc_out = "dfpt-phonon-o"
            self.files.tmp = "tmp"

            for i in range(self.ndtset+1):
                self.dataset[i].dfpt.status = True

            # overall default dataset
            self.dataset[0].electrons.params["getwfk"] = 1
            self.dataset[0].electrons.kpoints.params["kptopt"] = 3
            self.dataset[0].dfpt.params["rfphon"]  = 1
            self.dataset[0].dfpt.params["rfatpol"] = [1, self.dataset[0].system.xyz.natom]
            self.dataset[0].dfpt.params["rfdir"] = [1, 1, 1]
            self.dataset[0].electrons.use_tol(tol="tolvrs", value=1.0e-8)

            # 1) ground state scf calculation
            self.dataset[1].electrons.use_tol(tol="tolvrs", value=1.0e-10)
            self.dataset[1].electrons.kpoints.params["ngkpt"] = self.dataset[0].electrons.kpoints.params["ngkpt"]


            # 2) calculate ddk wave functions (iscf = -3, rfelfd =2, qpt = 0 0 0, rfdir = 1 1 1)
            self.dataset[2].electrons.use_tol(tol="tolwfr", value=1.0e-22)
            self.dataset[2].electrons.params["getwfk"] = -1
            self.dataset[2].electrons.kpoints.params["kptopt"] = 2
            self.dataset[2].electrons.kpoints.params["ngkpt"] = self.dataset[0].electrons.kpoints.params["ngkpt"]
            self.dataset[2].electrons.params["iscf"] = -3
            self.dataset[2].dfpt.params["rfelfd"] = 2
            self.dataset[2].dfpt.params["rfphon"] = 0
            self.dataset[2].dfpt.params["nqpt"] = 1
            self.dataset[2].dfpt.params["qpt"] = [0, 0, 0]


            # 3) Response function calculation of Q=0 phonons and electric field pert.

            self.dataset[3].electrons.params["getddk"] = 2
            self.dataset[3].electrons.kpoints.params["kptopt"] = 2
            self.dataset[3].electrons.kpoints.params["ngkpt"] = self.dataset[0].electrons.kpoints.params["ngkpt"]
            self.dataset[3].dfpt.params["rfelfd"] = 3
            self.dataset[3].dfpt.params["nqpt"] = 1
            self.dataset[3].dfpt.params["qpt"] = [0, 0, 0]
            self.dataset[3].electrons.use_tol(tol="tolvrs", value=1.0e-8)

            # calculate for q point other than gamma in qpath
            for i in range(4, self.ndtset+1):
                specialk = []
                for kpoint in self.qpath:
                    if kpoint[3] not in specialk:
                        self.dataset[i].dfpt.params["qpt"] = kpoint[0:3]
                        self.dataset[i].electrons.use_tol(tol="tolvrs", value=1.0e-8)
                        self.dataset[i].electrons.kpoints.params["ngkpt"] = self.dataset[0].electrons.kpoints.params["ngkpt"]
                        specialk.append(kpoint[3])
                    else:
                        pass
            #
            # generate pbs job submit script
            self.gen_pbs(directory=directory, script="dfpt-phonon.pbs", cmd="abinit", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # generate local bash job run script
            self.gen_bash(directory=directory, script="dfpt-phonon.sh", cmd="abinit", mpi=self.run_params["mpi"])

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash dfpt-phonon.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="dfpt-phonon", server=self.run_params["server"])
