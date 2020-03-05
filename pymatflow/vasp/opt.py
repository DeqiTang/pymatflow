#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil

from pymatflow.remote.server import server_handle
from pymatflow.vasp.vasp import vasp

"""
usage:
"""

class opt_run(vasp):
    """
    """
    def __init__(self):
        super().__init__()

        self.incar.set_runtype(runtype="opt")


    def optimize(self, directory="tmp-vasp-optimization", runopt="gen", auto=0):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.poscar.xyz.file, directory))

            with open(os.path.join(directory, "INCAR"), 'w') as fout:
                self.incar.to_incar(fout)
            with open(os.path.join(directory, "POSCAR"), 'w') as fout:
                self.poscar.to_poscar(fout)

            # gen slurm script
            self.gen_llhpc(directory=directory, cmd="vasp", scriptname="optimization.slurm")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="vasp_std", scriptname="optimization.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # gen local bash script
            self.gen_bash(directory=directory, cmd="vasp_std", scriptname="optimization.sh")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("vasp")
            os.system("bash optimization.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="optimization", server=self.run_params["server"])
    #


    def cubic(self, directory="tmp-vasp-opt-cubic", runopt="gen", auto=0, na=10, stepa=0.05):
        """
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
        os.system("cp %s %s/" % (self.poscar.xyz.file, directory))

        with open(os.path.join(directory, "KPOINTS"), "w") as fout:
            self.kpoints.to_kpoints(fout)

        with open(os.path.join(directory, "POSCAR"), 'w') as fout:
            self.poscar.to_poscar(fout)

        os.chdir(directory)

        # gen llhpc script
        with open("opt-cubic.slurm", 'w') as fout:
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

            a = self.poscar.xyz.cell[0][0]

            fout.write("v11=%f\n" % self.poscar.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.poscar.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.poscar.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.poscar.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.poscar.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.poscar.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.poscar.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.poscar.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.poscar.xyz.cell[2][2])

            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  mkdir relax-${a}\n")
            fout.write("  cp POTCAR KPOINTS INCAR relax-${a}/\n")
            fout.write("  cat > relax-${a}/POSCAR<<EOF\n")
            fout.write("general comment\n")
            fout.write("1.0\n")
            fout.write("${a} 0.000000 0.000000\n")
            fout.write("0.000000 ${a} 0.000000\n")
            fout.write("0.000000 0.000000 ${a}\n")
            fout.write("EOF\n")
            fout.write("  cat POSCAR | tail -n +6 >> relax-${a}/POSCAR\n")
            fout.write("  cd relax-${a}/\n")
            fout.write("  yhrun $PMF_VASP_STD\n")
            fout.write("  cd ../\n")
            fout.write("done\n")


        # gen pbs script
        with open("opt-cubic.pbs", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % self.run_params["jobname"])
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            fout.write("cat > INCAR<<EOF\n")
            self.incar.to_incar(fout)
            fout.write("EOF\n")

            a = self.poscar.xyz.cell[0][0]

            fout.write("v11=%f\n" % self.poscar.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.poscar.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.poscar.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.poscar.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.poscar.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.poscar.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.poscar.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.poscar.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.poscar.xyz.cell[2][2])

            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  mkdir relax-${a}\n")
            fout.write("  cp POTCAR KPOINTS INCAR relax-${a}/\n")
            fout.write("  cat > relax-${a}/POSCAR<<EOF\n")
            fout.write("general comment\n")
            fout.write("1.0\n")
            fout.write("${a} 0.000000 0.000000\n")
            fout.write("0.000000 ${a} 0.000000\n")
            fout.write("0.000000 0.000000 ${a}\n")
            fout.write("EOF\n")
            fout.write("  cat POSCAR | tail -n +6 >> relax-${a}/POSCAR\n")
            fout.write("  cd relax-${a}/\n")
            fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPMI_FABRICS shm:tmi vasp_std\n")
            fout.write("  cd ../\n")
            fout.write("done\n")


        # gen local bash script
        with open("opt-cubic.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")

            fout.write("cat > INCAR<<EOF\n")
            self.incar.to_incar(fout)
            fout.write("EOF\n")

            a = self.poscar.xyz.cell[0][0]

            fout.write("v11=%f\n" % self.poscar.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.poscar.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.poscar.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.poscar.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.poscar.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.poscar.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.poscar.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.poscar.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.poscar.xyz.cell[2][2])

            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  mkdir relax-${a}\n")
            fout.write("  cp POTCAR KPOINTS INCAR relax-${a}/\n")
            fout.write("  cat > relax-${a}/POSCAR<<EOF\n")
            fout.write("general comment\n")
            fout.write("1.0\n")
            fout.write("${a} 0.000000 0.000000\n")
            fout.write("0.000000 ${a} 0.000000\n")
            fout.write("0.000000 0.000000 ${a}\n")
            fout.write("EOF\n")
            fout.write("  cat POSCAR | tail -n +6 >> relax-${a}/POSCAR\n")
            fout.write("  cd relax-${a}/\n")
            fout.write("  %s $PMF_VASP_STD\n" % self.run_params["mpi"])
            fout.write("  cd ../\n")
            fout.write("done\n")

        # generate result analysis script
        os.system("mkdir -p post-processing")

        with open("post-processing/get_energy.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")
            # the comment
            fout.write("cat > energy-latconst.data <<EOF\n")
            fout.write("# format: a energy(eV)\n")
            fout.write("EOF\n")
            # end
            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  energy=`cat ../relax-${a}/OUTCAR | grep 'energy  without entropy=' | tail -1`\n")
            fout.write("  cat >> energy-latconst.data <<EOF\n")
            #fout.write("${a} ${energy:27:-36}\n")
            fout.write("${a} ${energy:27:17}\n")
            fout.write("EOF\n")
            fout.write("done\n")
            fout.write("cat > energy-latconst.gp<<EOF\n")
            fout.write("set term gif\n")
            fout.write("set output 'energy-latconst.gif'\n")
            fout.write("set title 'Energy Latconst'\n")
            fout.write("set xlabel 'latconst(a)'\n")
            fout.write("set ylabel 'Energy'\n")
            fout.write("plot 'energy-latconst.data' w l\n")
            fout.write("EOF\n")
            fout.write("\n")
            fout.write("gnuplot energy-latconst.gp")

        #os.system("cd post-processing; bash get_energy.sh; cd ../")
        os.chdir("../")
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash opt-cubic.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="opt-cubic", server=self.run_params["server"])


    def hexagonal(self, directory="tmp-vasp-opt-hexagonal", runopt="gen", auto=0, na=10, nc=10, stepa=0.05, stepc=0.05):
        """
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
        os.system("cp %s %s/" % (self.poscar.xyz.file, directory))

        with open(os.path.join(directory, "KPOINTS"), "w") as fout:
            self.kpoints.to_kpoints(fout)
        with open(os.path.join(directory, "POSCAR"), 'w') as fout:
            self.poscar.to_poscar(fout)

        os.chdir(directory)

        # gen llhpc script
        with open("opt-hexagonal.slurm", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
            fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
            fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
            fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
            fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
            fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])

            a = self.poscar.xyz.cell[0][0]
            c = self.poscar.xyz.cell[2][2]
            fout.write("v11=%f\n" % self.poscar.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.poscar.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.poscar.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.poscar.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.poscar.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.poscar.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.poscar.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.poscar.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.poscar.xyz.cell[2][2])
            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${a}-${c}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${a}-${c}/\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  cat > relax-${a}-${c}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${a}-${c}/POSCAR\n")
                    fout.write("  cd relax-${a}-${c}/\n")
                    fout.write("  yhrun $PMF_VASP_STD\n")
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  mkdir relax-${a}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${a}/\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  cat > relax-${a}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${a}/POSCAR\n")
                    fout.write("  cd relax-${a}/\n")
                    fout.write("  yhrun $PMF_VASP_STD\n")
                    fout.write("  cd ../\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${c}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${c}/\n")
                    fout.write("  cat > relax-${c}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("${v21} ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${c}/POSCAR\n")
                    fout.write("  cd relax-${c}/\n")
                    fout.write("  yhrun $PMF_VASP_STD\n")
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass


        # gen pbs script
        with open("opt-hexagonal.pbs", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % self.run_params["jobname"])
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("cat > INCAR<<EOF\n")
            self.incar.to_incar(fout)
            fout.write("EOF\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")

            a = self.poscar.xyz.cell[0][0]
            c = self.poscar.xyz.cell[2][2]
            fout.write("v11=%f\n" % self.poscar.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.poscar.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.poscar.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.poscar.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.poscar.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.poscar.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.poscar.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.poscar.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.poscar.xyz.cell[2][2])
            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${a}-${c}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${a}-${c}/\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  cat > relax-${a}-${c}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${a}-${c}/POSCAR\n")
                    fout.write("  cd relax-${a}-${c}/\n")
                    fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPMI_FABRICS shm:tmi $PMF_VASP_STD\n")
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  mkdir relax-${a}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${a}/\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  cat > relax-${a}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${a}/POSCAR\n")
                    fout.write("  cd relax-${a}/\n")
                    fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPMI_FABRICS shm:tmi $PMF_VASP_STD\n")
                    fout.write("  cd ../\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${c}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${c}/\n")
                    fout.write("  cat > relax-${c}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("${v21} ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${c}/POSCAR\n")
                    fout.write("  cd relax-${c}/\n")
                    fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPMI_FABRICS shm:tmi $PMF_VASP_STD\n")
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass

        # gen local bash script
        with open("opt-hexagonal.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")

            fout.write("cat > INCAR<<EOF\n")
            self.incar.to_incar(fout)
            fout.write("EOF\n")

            a = self.poscar.xyz.cell[0][0]
            c = self.poscar.xyz.cell[2][2]
            fout.write("v11=%f\n" % self.poscar.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.poscar.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.poscar.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.poscar.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.poscar.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.poscar.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.poscar.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.poscar.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.poscar.xyz.cell[2][2])
            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${a}-${c}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${a}-${c}/\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  cat > relax-${a}-${c}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${a}-${c}/POSCAR\n")
                    fout.write("  cd relax-${a}-${c}/\n")
                    fout.write("  %s $PMF_VASP_STD\n" % self.run_params["mpi"])
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  mkdir relax-${a}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${a}/\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  cat > relax-${a}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${a}/POSCAR\n")
                    fout.write("  cd relax-${a}/\n")
                    fout.write("  %s $PMF_VASP_STD\n" % self.run_params["mpi"])
                    fout.write("  cd ../\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${c}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${c}/\n")
                    fout.write("  cat > relax-${c}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("${v21} ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${c}/POSCAR\n")
                    fout.write("  cd relax-${c}/\n")
                    fout.write("  %s $PMF_VASP_STD\n" % self.run_params["mpi"])
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass


        # generate result analysis script
        os.system("mkdir -p post-processing")

        with open("post-processing/get_energy.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")
            # the comment
            if na >= 2 and nc >= 2:
                fout.write("cat > energy-latconst.data <<EOF\n")
                fout.write("# format: a c energy(eV)\n")
                fout.write("EOF\n")
            if na >= 2 and nc < 2:
                fout.write("cat > energy-latconst.data <<EOF\n")
                fout.write("# format: a energy(eV)\n")
                fout.write("EOF\n")
            if na < 2 and nc >= 2:
                fout.write("cat > energy-latconst.data <<EOF\n")
                fout.write("# format: c energy(eV)\n")
                fout.write("EOF\n")
            # end
            if na >= 2:
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  energy=`cat ../relax-${a}-${c}/OUTCAR | grep 'energy  without entropy=' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    #fout.write("${a} ${c} ${energy:27:-36}\n")
                    fout.write("${a} ${c} ${energy:27:17}\n")
                    fout.write("EOF\n")
                    fout.write("done\n")
                    fout.write("done\n")

                    fout.write("cat > energy-latconst.gp<<EOF\n")
                    fout.write("set term gif\n")
                    fout.write("set output 'energy-latconst.gif'\n")
                    fout.write("set title 'Energy Latconst'\n")
                    fout.write("set xlabel 'latconst(a)'\n")
                    fout.write("set ylabel 'latconst(c)'\n")
                    fout.write("set zlabel 'Energy'\n")
                    fout.write("splot 'energy-latconst.data'\n")
                    fout.write("EOF\n")
                    fout.write("\n")
                    fout.write("gnuplot energy-latconst.gp")
                else:
                    fout.write("  energy=`cat ../relax-${a}/OUTCAR | grep 'energy  without entropy=' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    #fout.write("${a} ${energy:27:-36}\n")
                    fout.write("${a} ${energy:27:17}\n")
                    fout.write("EOF\n")
                    fout.write("done\n")

                    fout.write("cat > energy-latconst.gp<<EOF\n")
                    fout.write("set term gif\n")
                    fout.write("set output 'energy-latconst.gif'\n")
                    fout.write("set title 'Energy Latconst'\n")
                    fout.write("set xlabel 'latconst(a)'\n")
                    fout.write("set ylabel 'Energy'\n")
                    fout.write("plot 'energy-latconst.data' w l\n")
                    fout.write("EOF\n")
                    fout.write("\n")
                    fout.write("gnuplot energy-latconst.gp")
            else:
                if nc >= 2:
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  energy=`cat ../relax-${c}/OUTCAR | grep 'energy  without entropy=' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    #fout.write("${c} ${energy:27:-36}\n")
                    fout.write("${c} ${energy:27:17}\n")
                    fout.write("EOF\n")
                    fout.write("done\n")
                    fout.write("cat > energy-latconst.gp<<EOF\n")
                    fout.write("set term gif\n")
                    fout.write("set output 'energy-latconst.gif'\n")
                    fout.write("set title 'Energy Latconst'\n")
                    fout.write("set xlabel 'latconst(c)'\n")
                    fout.write("set ylabel 'Energy'\n")
                    fout.write("plot 'energy-latconst.data' w l\n")
                    fout.write("EOF\n")
                    fout.write("\n")
                    fout.write("gnuplot energy-latconst.gp")
                else:
                    # nothing to do
                    pass
        os.chdir("../")
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash opt-hexagonal.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="opt-hexagonal", server=self.run_params["server"])

    def tetragonal(self, directory="tmp-vasp-opt-tetragonal", runopt="gen", auto=0, na=10, nc=10, stepa=0.05, stepc=0.05):
        """
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
        os.system("cp %s %s/" % (self.poscar.xyz.file, directory))


        with open(os.path.join(directory, "KPOINTS"), 'w') as fout:
            self.kpoints.to_kpoints(fout)
        with open(os.path.join(directory, "POSCAR"), 'w') as fout:
            self.poscar.to_poscar(fout)

        os.chdir(directory)

        # gen llhpc script
        with open("opt-tetragonal.slurm", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
            fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
            fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
            fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
            fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
            fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])

            a = self.poscar.xyz.cell[0][0]
            c = self.poscar.xyz.cell[2][2]

            fout.write("v11=%f\n" % self.poscar.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.poscar.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.poscar.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.poscar.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.poscar.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.poscar.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.poscar.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.poscar.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.poscar.xyz.cell[2][2])

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${a}-${c}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${a}-${c}/\n")
                    fout.write("  cat > relax-${a}-${c}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${a}-${c}/POSCAR\n")
                    fout.write("  cd relax-${a}-${c}/\n")
                    fout.write("  yhrun $PMF_VASP_STD\n")
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  mkdir relax-${a}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${a}/\n")
                    fout.write("  cat > relax-${a}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${a}/POSCAR\n")
                    fout.write("  cd relax-${a}/\n")
                    fout.write("  yhrun $PMF_VASP_STD\n")
                    fout.write("  cd ../\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${c}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${c}/\n")
                    fout.write("  cat > relax-${c}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("0.000000 ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${c}/POSCAR\n")
                    fout.write("  cd relax-${c}/\n")
                    fout.write("  yhrun $PMF_VASP_STD\n")
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass


        # gen pbs script
        with open("opt-tetragonal.pbs", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % self.run_params["jobname"])
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("cat > INCAR<<EOF\n")
            self.incar.to_incar(fout)
            fout.write("EOF\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")

            a = self.poscar.xyz.cell[0][0]
            c = self.poscar.xyz.cell[2][2]

            fout.write("v11=%f\n" % self.poscar.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.poscar.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.poscar.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.poscar.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.poscar.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.poscar.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.poscar.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.poscar.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.poscar.xyz.cell[2][2])

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${a}-${c}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${a}-${c}/\n")
                    fout.write("  cat > relax-${a}-${c}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${a}-${c}/POSCAR\n")
                    fout.write("  cd relax-${a}-${c}/\n")
                    fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPMI_FABRICS shm:tmi $PMF_VASP_STD\n")
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  mkdir relax-${a}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${a}/\n")
                    fout.write("  cat > relax-${a}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${a}/POSCAR\n")
                    fout.write("  cd relax-${a}/\n")
                    fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPMI_FABRICS shm:tmi $PMF_VASP_STD\n")
                    fout.write("  cd ../\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${c}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${c}/\n")
                    fout.write("  cat > relax-${c}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("0.000000 ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${c}/POSCAR\n")
                    fout.write("  cd relax-${c}/\n")
                    fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPMI_FABRICS shm:tmi $PMF_VASP_STD\n")
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass

        # gen local bash script
        with open("opt-tetragonal.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")

            fout.write("cat > INCAR<<EOF\n")
            self.incar.to_incar(fout)
            fout.write("EOF\n")

            a = self.poscar.xyz.cell[0][0]
            c = self.poscar.xyz.cell[2][2]

            fout.write("v11=%f\n" % self.poscar.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.poscar.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.poscar.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.poscar.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.poscar.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.poscar.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.poscar.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.poscar.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.poscar.xyz.cell[2][2])

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${a}-${c}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${a}-${c}/\n")
                    fout.write("  cat > relax-${a}-${c}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${a}-${c}/POSCAR\n")
                    fout.write("  cd relax-${a}-${c}/\n")
                    fout.write("  %s $PMF_VASP_STD\n" % self.run_params["mpi"])
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  mkdir relax-${a}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${a}/\n")
                    fout.write("  cat > relax-${a}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${a}/POSCAR\n")
                    fout.write("  cd relax-${a}/\n")
                    fout.write("  %s $PMF_VASP_STD\n" % self.run_params["mpi"])
                    fout.write("  cd ../\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${c}\n")
                    fout.write("  cp POTCAR KPOINTS INCAR relax-${c}/\n")
                    fout.write("  cat > relax-${c}/POSCAR<<EOF\n")
                    fout.write("general comment\n")
                    fout.write("1.0\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("0.000000 ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat POSCAR | tail -n +6 >> relax-${c}/POSCAR\n")
                    fout.write("  cd relax-${c}/\n")
                    fout.write("  %s $PMF_VASP_STD\n" % self.run_params["mpi"])
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass

        # generate result analysis script
        os.system("mkdir -p post-processing")

        with open("post-processing/get_energy.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")
            # the comment
            if na >= 2 and nc >= 2:
                fout.write("cat > energy-latconst.data <<EOF\n")
                fout.write("# format: a c energy(eV)\n")
                fout.write("EOF\n")
            if na >= 2 and nc < 2:
                fout.write("cat > energy-latconst.data <<EOF\n")
                fout.write("# format: a energy(eV)\n")
                fout.write("EOF\n")
            if na < 2 and nc >= 2:
                fout.write("cat > energy-latconst.data <<EOF\n")
                fout.write("# format: c energy(eV)\n")
                fout.write("EOF\n")
            # end
            if na >= 2:
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  energy=`cat ../relax-${a}-${c}/OUTCAR | grep 'energy  without entropy=' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    #fout.write("${a} ${c} ${energy:27:-36}\n")
                    fout.write("${a} ${c} ${energy:27:17}\n")
                    fout.write("EOF\n")
                    fout.write("done\n")
                    fout.write("done\n")

                    fout.write("cat > energy-latconst.gp<<EOF\n")
                    fout.write("set term gif\n")
                    fout.write("set output 'energy-latconst.gif'\n")
                    fout.write("set title 'Energy Latconst'\n")
                    fout.write("set xlabel 'latconst(a)'\n")
                    fout.write("set ylabel 'latconst(c)'\n")
                    fout.write("set zlabel 'Energy'\n")
                    fout.write("splot 'energy-latconst.data'\n")
                    fout.write("EOF\n")
                    fout.write("\n")
                    fout.write("gnuplot energy-latconst.gp")
                else:
                    fout.write("  energy=`cat ../relax-${a}/OUTCAR | grep 'energy  without entropy=' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    #fout.write("${a} ${energy:27:-36}\n")
                    fout.write("${a} ${energy:27:17}\n")
                    fout.write("EOF\n")
                    fout.write("done\n")

                    fout.write("cat > energy-latconst.gp<<EOF\n")
                    fout.write("set term gif\n")
                    fout.write("set output 'energy-latconst.gif'\n")
                    fout.write("set title 'Energy Latconst'\n")
                    fout.write("set xlabel 'latconst(a)'\n")
                    fout.write("set ylabel 'Energy'\n")
                    fout.write("plot 'energy-latconst.data' w l\n")
                    fout.write("EOF\n")
                    fout.write("\n")
                    fout.write("gnuplot energy-latconst.gp")
            else:
                if nc >= 2:
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  energy=`cat ../relax-${c}/OUTCAR | grep 'energy  without entropy=' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    #fout.write("${c} ${energy:27:-36}\n")
                    fout.write("${c} ${energy:27:17}\n")
                    fout.write("EOF\n")
                    fout.write("done\n")
                    fout.write("cat > energy-latconst.gp<<EOF\n")
                    fout.write("set term gif\n")
                    fout.write("set output 'energy-latconst.gif'\n")
                    fout.write("set title 'Energy Latconst'\n")
                    fout.write("set xlabel 'latconst(c)'\n")
                    fout.write("set ylabel 'Energy'\n")
                    fout.write("plot 'energy-latconst.data' w l\n")
                    fout.write("EOF\n")
                    fout.write("\n")
                    fout.write("gnuplot energy-latconst.gp")
                else:
                    # nothing to do
                    pass
        #os.system("cd post-processing; bash get_energy.sh; cd ../")
        os.chdir("../")
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash opt-tetragonal.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="opt-tetragonal", server=self.run_params["server"])
