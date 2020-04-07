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

Note:
    Actually there are many combination of response function calculation implemented in Abinit.
    For instance, there are ways to  calculate Born effective charge. you can calculate it
    from phonon response or from electric field response.
    we know rfelfd=2 we can do ddk calculation, and ddk can also be calculated using brrryopt=-2

"""
import os
import sys
import shutil


from pymatflow.remote.server import server_handle
from pymatflow.abinit.abinit import abinit

from pymatflow.abinit.base.guard import abinit_guard


class dfpt_elastic_piezo_dielec(abinit):
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

    def run(self, directory="tmp-abinit-dfpt-elastic-piezo-dielec", runopt="gen", auto=0):
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
            self.dataset[2].dfpt.params["rfelfd"] = 2 # here we only do ddk calculation, no electric field perturbation
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
            # in dataset 2 we only do the ddk calculation, and int the output file you can only find Effective charge
            # (from phonon response). now we can do electric field perturbation using rfelfd=3,
            # this way Born effective charge (from electric fiedl perturbation) will also be print out to the abinit output
            # and output of anaddb will contains information of Effective charge.
            # plus the output can be analysed further by anaddb and dielectric tensor can be calculated.
            self.dataset[3].dfpt.params["rfelfd"] = 3 # you can also set it to None, if you want no perturbation to electric field
            self.dataset[3].dfpt.params["rfatpol"] = [1, self.dataset[0].system.xyz.natom]
            self.dataset[3].dfpt.params["rfdir"] = [1, 1, 1]
            self.dataset[3].dfpt.params["rfstrs"] = 3

            # llhpc jobsubmit script
            with open(os.path.join(directory, "dfpt-elastic-piezo-dielec.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                fout.write("cat > %s<<EOF\n" % self.files.main_in)
                fout.write(self.to_string())
                fout.write("EOF\n")
                fout.write("cat > %s<<EOF\n" % self.files.name)
                fout.write(self.files.to_string(system=self.dataset[0].system))
                fout.write("EOF\n")
                fout.write("yhrun %s < %s\n" % ("$PMF_ABINIT", self.files.name))

                # use anaddb to analyse DDB file generated previously
                #with open(os.path.join(directory, "anaddb.in"), "w") as fout:
                fout.write("cat > %s <<EOF\n" % "anaddb.in")
                fout.write("elaflag 3\n")
                fout.write("piezoflag 7\n")
                # piezoflag = 7 means calculate all the possible piezoelectric tensors, including e (clamped and relaxed ion), d, g and h tensors
                # if g and h tensor are to be calculated must set dieflag to 3 or 4, or it will print the wrong value
                # but the calculation will continue
                fout.write("instrflag 1\n")
                fout.write("chneut 1\n")
                fout.write("dieflag 1\n")
                fout.write("# if ecut and kmesh is not sufficient, the calc of dielectric tensor might collapse\n")
                #  The frequency-dependent dielectric tensor is calculated
                # see https://docs.abinit.org/variables/anaddb/#dieflag
                # if ecut and kmesh is not sufficient, the calc of dielectric tensor might collapse
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
                fout.write("yhrun %s < %s\n" % ("$PMF_ANADDB", "anaddb.files"))


            # pbs jobsubmit script
            with open(os.path.join(directory, "dfpt-elastic-piezo-dielec.pbs"), 'w') as fout:
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
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s\n" % ("$PMF_ABINIT", self.files.name))

                # use anaddb to analyse DDB file generated previously
                #with open(os.path.join(directory, "anaddb.in"), "w") as fout:
                fout.write("cat > %s <<EOF\n" % "anaddb.in")
                fout.write("elaflag 3\n")
                fout.write("piezoflag 7\n")
                # piezoflag = 7 means calculate all the possible piezoelectric tensors, including e (clamped and relaxed ion), d, g and h tensors
                # if g and h tensor are to be calculated must set dieflag to 3 or 4, or it will print the wrong value
                # but the calculation will continue
                fout.write("instrflag 1\n")
                fout.write("chneut 1\n")
                fout.write("dieflag 1\n")
                fout.write("# if ecut and kmesh is not sufficient, the calc of dielectric tensor might collapse\n")
                #  The frequency-dependent dielectric tensor is calculated
                # see https://docs.abinit.org/variables/anaddb/#dieflag
                # if ecut and kmesh is not sufficient, the calc of dielectric tensor might collapse
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
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s\n" % ("$PMF_ANADDB", "anaddb.files"))

            # local bash script
            with open(os.path.join(directory, "dfpt-elastic-piezo-dielec.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n")

                fout.write("cat > %s<<EOF\n" % self.files.main_in)
                fout.write(self.to_string())
                fout.write("EOF\n")
                fout.write("cat > %s<<EOF\n" % self.files.name)
                fout.write(self.files.to_string(system=self.dataset[0].system))
                fout.write("EOF\n")
                fout.write("%s %s < %s\n" % (self.run_params["mpi"], "$PMF_ABINIT", self.files.name))

                # use anaddb to analyse DDB file generated previously
                #with open(os.path.join(directory, "anaddb.in"), "w") as fout:
                fout.write("cat > %s <<EOF\n" % "anaddb.in")
                fout.write("elaflag 3\n")
                fout.write("piezoflag 7 \n")
                # piezoflag = 7 means calculate all the possible piezoelectric tensors, including e (clamped and relaxed ion), d, g and h tensors
                # if g and h tensor are to be calculated must set dieflag to 3 or 4, or it will print the wrong value
                # but the calculation will continue
                fout.write("instrflag 1\n")
                fout.write("chneut 1\n")
                fout.write("dieflag 1\n")
                fout.write("# if ecut and kmesh is not sufficient, the calc of dielectric tensor might collapse\n")
                #  The frequency-dependent dielectric tensor is calculated
                # see https://docs.abinit.org/variables/anaddb/#dieflag
                # if ecut and kmesh is not sufficient, the calc of dielectric tensor might collapse
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
                fout.write("%s %s < %s\n" % (self.run_params["mpi"], "$PMF_ANADDB", "anaddb.files"))

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash dfpt-elastic-piezo-dielec.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="dfpt-elastic-piezo-dielec", server=self.run_params["server"])




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

            self.dataset[0].electrons.kpoints.params["kptopt"] = 3
            #self.dataset[0].dfpt.params["rfphon"]  = 1
            # If one of rfphon, rfddk, rfelfd, or rfstrs is non-zero, while optdriver is not defined in the input file,
            # ABINIT will set optdriver to 1 automatically(means for response function calc)
            # so we should not set them in the dataset 0 for the default setting.
            # but rfatpol and rfdir can be set in the dataset 0 which will not affect default value of optdriver
            self.dataset[0].dfpt.params["rfatpol"] = [1, self.dataset[0].system.xyz.natom]
            self.dataset[0].dfpt.params["rfdir"] = [1, 1, 1]
            self.dataset[0].electrons.use_tol(tol="tolvrs", value=1.0e-8)

            # 1) ground state scf calculation
            self.dataset[1].electrons.params["iscf"] = 7
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
            self.dataset[3].electrons.params["getwfk"] = 1
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
                        self.dataset[i].electrons.params["getwfk"] = 1
                        self.dataset[i].dfpt.params["qpt"] = kpoint[0:3]
                        self.dataset[i].electrons.use_tol(tol="tolvrs", value=1.0e-8)
                        self.dataset[i].electrons.kpoints.params["ngkpt"] = self.dataset[0].electrons.kpoints.params["ngkpt"]
                        specialk.append(kpoint[3])
                    else:
                        pass
            #
            # generate pbs job submit script
            #self.gen_pbs(directory=directory, script="dfpt-phonon.pbs", cmd="$PMF_ABINIT", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # generate local bash job run script
            #self.gen_bash(directory=directory, script="dfpt-phonon.sh", cmd="$PMF_ABINIT", mpi=self.run_params["mpi"])

            # llhpc jobsubmit script
            with open(os.path.join(directory, "dfpt-phonon.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                fout.write("cat > %s<<EOF\n" % self.files.main_in)
                fout.write(self.to_string())
                fout.write("EOF\n")
                fout.write("cat > %s<<EOF\n" % self.files.name)
                fout.write(self.files.to_string(system=self.dataset[0].system))
                fout.write("EOF\n")
                fout.write("yhrun %s < %s\n" % ("$PMF_ABINIT", self.files.name))

                # use mrgddb and anaddb to analyse DDB file generated previously
                fout.write("cat %s <<EOF\n" % "mrgddb.in")
                fout.write("mrgddb.ddb.out\n")
                fout.write("xxx\n")
                fout.write("%d\n" % (self.ndtset-2))
                for i in range(3, self.ndtset+1):
                    fout.write("dfpt-phonon-o_DS%d_DDB\n" % i)
                fout.write("EOF\n")
                fout.write("yhrun $PMF_MRGDDB < mrgddb.in\n")

                # anaddb
                fout.write("cat > %s <<EOF\n" % "anaddb.in")
                fout.write("ifcflag 1\n")
                fout.write("ifcout 0\n")
                fout.write("! wavevector grid\n")
                fout.write("brav 2\n")
                fout.write("ngqpt 3 3 3\n")
                fout.write("! effective charge\n")
                fout.write("chneut 1\n")
                fout.write("! interatomic force constant info\n")
                fout.write("dipdip 1\n")
                fout.write("! Phonon band structure output for band2eps\n")
                fout.write("eivec 4\n")
                fout.write("! wavevector list\n")



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
                #fout.write("yhrun %s < %s\n" % ("$PMF_ANADDB", "anaddb.files"))


            # pbs jobsubmit script
            with open(os.path.join(directory, "dfpt-phonon.pbs"), 'w') as fout:
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
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s\n" % ("$PMF_ABINIT", self.files.name))

                # use mrgddb and anaddb to analyse DDB file generated previously
                fout.write("cat %s <<EOF\n" % "mrgddb.in")
                fout.write("mrgddb.ddb.out\n")
                fout.write("xxx\n")
                fout.write("%d\n" % (self.ndtset-2))
                for i in range(3, self.ndtset+1):
                    fout.write("dfpt-phonon-o_DS%d_DDB\n" % i)
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_MRGDDB < mrgddb.in\n")

                # anaddb
                fout.write("cat > %s <<EOF\n" % "anaddb.in")
                fout.write("ifcflag 1\n")
                fout.write("ifcout 0\n")
                fout.write("! wavevector grid\n")
                fout.write("brav 2\n")
                fout.write("ngqpt 3 3 3\n")
                fout.write("! effective charge\n")
                fout.write("chneut 1\n")
                fout.write("! interatomic force constant info\n")
                fout.write("dipdip 1\n")
                fout.write("! Phonon band structure output for band2eps\n")
                fout.write("eivec 4\n")
                fout.write("! wavevector list\n")



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
                #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s\n" % ("$PMF_ANADDB", "anaddb.files"))

            # local bash script
            with open(os.path.join(directory, "dfpt-phonon.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n")

                fout.write("cat > %s<<EOF\n" % self.files.main_in)
                fout.write(self.to_string())
                fout.write("EOF\n")
                fout.write("cat > %s<<EOF\n" % self.files.name)
                fout.write(self.files.to_string(system=self.dataset[0].system))
                fout.write("EOF\n")
                fout.write("%s %s < %s\n" % (self.run_params["mpi"], "$PMF_ABINIT", self.files.name))

                # use mrgddb and anaddb to analyse DDB file generated previously
                fout.write("cat %s <<EOF\n" % "mrgddb.in")
                fout.write("mrgddb.ddb.out\n")
                fout.write("xxx\n")
                fout.write("%d\n" % (self.ndtset-2))
                for i in range(3, self.ndtset+1):
                    fout.write("dfpt-phonon-o_DS%d_DDB\n" % i)
                fout.write("EOF\n")
                fout.write("%s $PMF_MRGDDB < mrgddb.in\n" % self.run_params["mpi"])

                # anaddb
                fout.write("cat > %s <<EOF\n" % "anaddb.in")
                fout.write("ifcflag 1\n")
                fout.write("ifcout 0\n")
                fout.write("! wavevector grid\n")
                fout.write("brav 2\n")
                fout.write("ngqpt 3 3 3\n")
                fout.write("! effective charge\n")
                fout.write("chneut 1\n")
                fout.write("! interatomic force constant info\n")
                fout.write("dipdip 1\n")
                fout.write("! Phonon band structure output for band2eps\n")
                fout.write("eivec 4\n")
                fout.write("! wavevector list\n")



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
                #fout.write("%s %s < %s\n" % (self.run_params["mpi"], "$PMF_ANADDB", "anaddb.files"))


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash dfpt-phonon.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="dfpt-phonon", server=self.run_params["server"])
