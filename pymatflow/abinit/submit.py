import os 

from pymatflow.remote.server import server_handle


def submit(abinit, directory, prefix, runopt, auto, optic=None):
    # llhpc jobsubmit script
    with open(os.path.join(directory, "%s.slurm" % prefix), 'w') as fout:
        fout.write("#!/bin/bash\n")
        fout.write("#SBATCH -p %s\n" % abinit.run_params["partition"])
        fout.write("#SBATCH -N %d\n" % abinit.run_params["nodes"])
        fout.write("#SBATCH -n %d\n" % abinit.run_params["ntask"])
        fout.write("#SBATCH -J %s\n" % abinit.run_params["jobname"])
        fout.write("#SBATCH -o %s\n" % abinit.run_params["stdout"])
        fout.write("#SBATCH -e %s\n" % abinit.run_params["stderr"])
        #fout.write("cat > %s<<EOF\n" % abinit.files.main_in)
        fout.write("cat > %s.abi<<EOF\n" % prefix)
        fout.write(abinit.to_string())
        fout.write("EOF\n")
        #fout.write("cat > %s<<EOF\n" % abinit.files.name)
        #fout.write(abinit.files.to_string(system=abinit.dataset[0].system))
        #fout.write("EOF\n")
        #fout.write("yhrun %s < %s\n" % ("$PMF_ABINIT", abinit.files.name))
        fout.write("yhrun %s %s.abi\n" % ("$PMF_ABINIT", prefix))

        if optic != None:
            fout.write("cat > %s <<EOF\n" % "optic.abi")
            fout.write(optic.to_string())
            fout.write("EOF\n")
            fout.write("%s %s\n" % ("$PMF_ABINIT_OPTIC", "optic.abi"))


    # pbs jobsubmit script
    with open(os.path.join(directory, "%s.pbs" % prefix), 'w') as fout:
        fout.write("#!/bin/bash\n")
        fout.write("#PBS -N %s\n" % abinit.run_params["jobname"])
        fout.write("#PBS -l nodes=%d:ppn=%d\n" % (abinit.run_params["nodes"], abinit.run_params["ppn"]))
        fout.write("\n")
        fout.write("cd $PBS_O_WORKDIR\n")
        fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
        #fout.write("cat > %s<<EOF\n" % abinit.files.main_in)
        fout.write("cat > %s.abi<<EOF\n" % prefix)
        fout.write(abinit.to_string())
        fout.write("EOF\n")
        #fout.write("cat > %s<<EOF\n" % abinit.files.name)
        #fout.write(abinit.files.to_string(system=abinit.dataset[0].system))
        #fout.write("EOF\n")
        #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s\n" % ("$PMF_ABINIT", abinit.files.name))
        fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s %s.abi\n" % ("$PMF_ABINIT", prefix))

        if optic != None:
            fout.write("cat > %s <<EOF\n" % "optic.abi")
            fout.write(optic.to_string())
            fout.write("EOF\n")
            fout.write("%s %s\n" % ("$PMF_ABINIT_OPTIC", "optic.abi"))


    # local bash script
    with open(os.path.join(directory, "%s.sh" % prefix), 'w') as fout:
        fout.write("#!/bin/bash\n")

        fout.write("cat > %s<<EOF\n" % abinit.files.main_in)
        fout.write(abinit.to_string())
        fout.write("EOF\n")
        #fout.write("cat > %s<<EOF\n" % abinit.files.name)
        fout.write("cat > %s.abi<<EOF\n" % prefix)
        fout.write(abinit.files.to_string(system=abinit.dataset[0].system))
        fout.write("EOF\n")
        #fout.write("%s %s < %s\n" % (abinit.run_params["mpi"], "$PMF_ABINIT", abinit.files.name))
        fout.write("%s %s %s.abi\n" % (abinit.run_params["mpi"], "$PMF_ABINIT", prefix))

        if optic != None:
            fout.write("cat > %s <<EOF\n" % "optic.abi")
            fout.write(optic.to_string())
            fout.write("EOF\n")
            fout.write("%s %s\n" % ("$PMF_ABINIT_OPTIC", "optic.abi"))

    # cdcloud jobsubmit script
    with open(os.path.join(directory, "%s.slurm_cd" % prefix), 'w') as fout:
        fout.write("#!/bin/bash\n")
        fout.write("#SBATCH -p %s\n" % abinit.run_params["partition"])
        fout.write("#SBATCH -N %d\n" % abinit.run_params["nodes"])
        fout.write("#SBATCH -n %d\n" % abinit.run_params["ntask"])
        fout.write("#SBATCH -J %s\n" % abinit.run_params["jobname"])
        fout.write("#SBATCH -o %s\n" % abinit.run_params["stdout"])
        fout.write("#SBATCH -e %s\n" % abinit.run_params["stderr"])
        fout.write("#\n")
        fout.write("export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so\n")
        fout.write("export FORT_BUFFERED=1\n")                
        #fout.write("cat > %s<<EOF\n" % abinit.files.main_in)
        fout.write("cat > %s.abi<<EOF\n" % prefix)
        fout.write(abinit.to_string())
        fout.write("EOF\n")
        #fout.write("cat > %s<<EOF\n" % abinit.files.name)
        #fout.write(abinit.files.to_string(system=abinit.dataset[0].system))
        #fout.write("EOF\n")
        #fout.write("srun --mpi=pmix_v3 %s < %s\n" % ("$PMF_ABINIT", abinit.files.name))
        fout.write("srun --mpi=pmix_v3 %s %s.abi\n" % ("$PMF_ABINIT", prefix))

        if optic != None:
            fout.write("cat > %s <<EOF\n" % "optic.abi")
            fout.write(optic.to_string())
            fout.write("EOF\n")
            fout.write("%s %s\n" % ("$PMF_ABINIT_OPTIC", "optic.abi"))


    if runopt == "run" or runopt == "genrun":
        os.chdir(directory)
        os.system("bash %s.sh" % prefix)
        os.chdir("../")
    server_handle(auto=auto, directory=directory, jobfilebase=prefix, server=abinit.run_params["server"])
