#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil


from pymatflow.remote.server import server_handle
from pymatflow.abinit.abinit import abinit

class opt_run(abinit):
    """
    """
    def __init__(self):
        super().__init__()
        self.dataset[0].electrons.basic_setting()
        self.dataset[0].ions.basic_setting(mode="opt")

        self.dataset[0].guard.set_queen(queen="opt")


    def optimize(self, directory="tmp-abinit-opt", runopt="gen", auto=0):

        self.dataset[0].electrons.set_scf_nscf("scf")

        self.files.name = "optimization.files"
        self.files.main_in = "optimization.in"
        self.files.main_out = "optimization.out"
        self.files.wavefunc_in = "optimization-i"
        self.files.wavefunc_out = "optimization-o"
        self.files.tmp = "tmp"

        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
            os.system("cp %s %s/" % (self.dataset[0].system.xyz.file, directory))

            #
            # generate llhpc job submit script
            self.gen_llhpc(directory=directory, script="optimization.slurm", cmd="$PMF_ABINIT")

            # generate pbs job submit script
            self.gen_pbs(directory=directory, script="optimization.pbs", cmd="$PMF_ABINIT", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # generate local bash job run script
            self.gen_bash(directory=directory, script="optimization.sh", cmd="$PMF_ABINIT", mpi=self.run_params["mpi"])


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.system("bash %s" % "optimization.sh")
            os.chdir("../")

        server_handle(auto=auto, directory=directory, jobfilebase="optimization", server=self.run_params["server"])

    def run(self, directory="tmp-abinit-opt", runopt="gen", auto=0):

        self.dataset[0].electrons.set_scf_nscf("scf")

        self.files.name = "optimization.files"
        self.files.main_in = "optimization.in"
        self.files.main_out = "optimization.out"
        self.files.wavefunc_in = "optimization-i"
        self.files.wavefunc_out = "optimization-o"
        self.files.tmp = "tmp"

        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
            os.system("cp %s %s/" % (self.dataset[0].system.xyz.file, directory))

            #
            # generate llhpc job submit script
            self.gen_llhpc(directory=directory, script="optimization.slurm", cmd="$PMF_ABINIT")
            # generate pbs job submit script
            self.gen_pbs(directory=directory, script="optimization.pbs", cmd="$PMF_ABINIT", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # generate local bash job run script
            self.gen_bash(directory=directory, script="optimization.sh", cmd="$PMF_ABINIT", mpi=self.run_params["mpi"])


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.system("bash %s" % "optimization.sh")
            os.chdir("../")

        server_handle(auto=auto, directory=directory, jobfilebase="optimization", server=self.run_params["server"])


    def cubic(self, directory="tmp-abinit-opt-cubic", runopt="genrun", auto=0, na=10, stepa=0.05):
        """
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        os.system("cp *.psp8 %s/" % directory)
        os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
        os.system("cp %s %s/" % (self.dataset[0].system.xyz.file, directory))

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
            fout.write("cat > optimization.in<<EOF\n")
            #self.dataset[0].to_input(fout)
            fout.write(self.dataset[0].to_string())
            fout.write("EOF\n")
            fout.write("cat > optimization.files<<EOF\n")
            #self.files.name = "optimization.files"
            self.files.main_in = "optimization.in"
            self.files.main_out = "optimization.out"
            self.files.wavefunc_in = "optimization-i"
            self.files.wavefunc_out = "optimization-o"
            self.files.tmp = "tmp"
            #self.files.to_files(fout, self.dataset[0].system)
            fout.write(self.files.to_string(system=self.dataset[0].system))
            fout.write("EOF\n")


            a = self.dataset[0].system.xyz.cell[0][0]

            fout.write("v11=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.dataset[0].system.xyz.cell[2][2])

            fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
            fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")
            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  mkdir relax-${a}\n")
            fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}/\n")
            fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}/optimization.in\n")
            fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
            fout.write("${a} 0.000000 0.000000\n")
            fout.write("0.000000 ${a} 0.000000\n")
            fout.write("0.000000 0.000000 ${a}\n")
            fout.write("EOF\n")
            fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}/optimization.in\n")
            fout.write("  cd relax-${a}/\n")
            fout.write("  yhrun $PMF_ABINIT < optimization.files\n")
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
            fout.write("cat > optimization.in<<EOF\n")
            #self.dataset[0].to_input(fout)
            fout.write(self.dataset[0].to_string())
            fout.write("EOF\n")
            fout.write("cat > optimization.files<<EOF\n")
            #self.files.name = "optimization.files"
            self.files.main_in = "optimization.in"
            self.files.main_out = "optimization.out"
            self.files.wavefunc_in = "optimization-i"
            self.files.wavefunc_out = "optimization-o"
            self.files.tmp = "tmp"
            #self.files.to_files(fout, self.dataset[0].system)
            fout.write(self.files.to_string(system=self.dataset[0].system))
            fout.write("EOF\n")


            a = self.dataset[0].system.xyz.cell[0][0]

            fout.write("v11=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.dataset[0].system.xyz.cell[2][2])

            fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
            fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")
            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  mkdir relax-${a}\n")
            fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}/\n")
            fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}/optimization.in\n")
            fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
            fout.write("${a} 0.000000 0.000000\n")
            fout.write("0.000000 ${a} 0.000000\n")
            fout.write("0.000000 0.000000 ${a}\n")
            fout.write("EOF\n")
            fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}/optimization.in\n")
            fout.write("  cd relax-${a}/\n")
            fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_ABINIT < optimization.files\n")
            fout.write("  cd ../\n")
            fout.write("done\n")

        # gen local bash script
        with open("opt-cubic.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")

            fout.write("cat > optimization.in<<EOF\n")
            #self.dataset[0].to_input(fout)
            fout.write(self.dataset[0].to_string())
            fout.write("EOF\n")
            fout.write("cat > optimization.files<<EOF\n")
            #self.files.name = "optimization.files"
            self.files.main_in = "optimization.in"
            self.files.main_out = "optimization.out"
            self.files.wavefunc_in = "optimization-i"
            self.files.wavefunc_out = "optimization-o"
            self.files.tmp = "tmp"
            #self.files.to_files(fout, self.dataset[0].system)
            fout.write(self.files.to_string(system=self.dataset[0].system))
            fout.write("EOF\n")


            a = self.dataset[0].system.xyz.cell[0][0]

            fout.write("v11=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.dataset[0].system.xyz.cell[2][2])

            fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
            fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")
            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  mkdir relax-${a}\n")
            fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}/\n")
            fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}/optimization.in\n")
            fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
            fout.write("${a} 0.000000 0.000000\n")
            fout.write("0.000000 ${a} 0.000000\n")
            fout.write("0.000000 0.000000 ${a}\n")
            fout.write("EOF\n")
            fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}/optimization.in\n")
            fout.write("  cd relax-${a}/\n")
            fout.write("  %s $PMF_ABINIT < optimization.files\n" % self.run_params["mpi"])
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
            fout.write("  energy=`cat ../relax-${a}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 1`\n" % self.files.main_out)
            fout.write("  cat >> energy-latconst.data <<EOF\n")
            fout.write("${a} ${energy}\n")
            fout.write("EOF\n")
            fout.write("done\n")
            fout.write("cat > energy-latconst.gp<<EOF\n")
            fout.write("set term gif\n")
            fout.write("set output 'energy-latconst.gif'\n")
            fout.write("set title Energy Latconst\n")
            fout.write("set xlabel 'latconst(a)'\n")
            fout.write("set ylabel 'Energy'\n")
            fout.write("plot 'energy-latconst.data' w l\n")
            fout.write("EOF\n")

        #os.system("cd post-processing; bash get_energy.sh; cd ../")
        os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash %s" % "opt-cubic.sh")
            os.chdir("../")

        server_handle(auto=auto, directory=directory, jobfilebase="opt-cubic", server=self.run_params["server"])

    def hexagonal(self, directory="tmp-abinit-opt-hexagonal", runopt="gen", auto=0, na=10, stepa=0.05, nc=10, stepc=0.05):
        """
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        os.system("cp *.psp8 %s/" % directory)
        os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
        os.system("cp %s %s/" % (self.dataset[0].system.xyz.file, directory))

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
            fout.write("cat > optimization.in<<EOF\n")
            #self.dataset[0].to_input(fout)
            fout.write(self.dataset[0].to_string())
            fout.write("EOF\n")
            fout.write("cat > optimization.files<<EOF\n")
            #self.files.name = "optimization.files"
            self.files.main_in = "optimization.in"
            self.files.main_out = "optimization.out"
            self.files.wavefunc_in = "optimization-i"
            self.files.wavefunc_out = "optimization-o"
            self.files.tmp = "tmp"
            #self.files.to_files(fout, self.dataset[0].system)
            fout.write(self.files.to_string(system=self.dataset[0].system))
            fout.write("EOF\n")

            a = self.dataset[0].system.xyz.cell[0][0]

            fout.write("v11=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.dataset[0].system.xyz.cell[2][2])

            fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
            fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")
            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${a}-${c}\n")
                    fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}-${c}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}-${c}/optimization.in\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  cat >> relax-${a}-${c}/optimization.in<<EOF\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}-${c}/optimization.in\n")
                    fout.write("  cd relax-${a}-${c}/\n")
                    fout.write("  yhrun $PMF_ABINIT < optimization.files\n")
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  mkdir relax-${a}\n")
                    fout.write("  cp optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}/optimization.in\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}/optimization.in\n")
                    fout.write("  cd relax-${a}/\n")
                    fout.write("  yhrun $PMF_ABINIT < optimization.files\n")
                    fout.write("  cd ../\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${c}\n")
                    fout.write("  cp optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${c}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${c}/optimization.in\n")
                    fout.write("  cat >> relax-${c}/optimization.in<<EOF\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("${v21} ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${c}/optimization.in\n")
                    fout.write("  cd relax-${c}/\n")
                    fout.write("  yhrun $PMF_ABINIT < optimization.files\n")
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
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            fout.write("cat > optimization.in<<EOF\n")
            #self.dataset[0].to_input(fout)
            fout.write(self.dataset[0].to_string())
            fout.write("EOF\n")
            fout.write("cat > optimization.files<<EOF\n")
            #self.files.name = "optimization.files"
            self.files.main_in = "optimization.in"
            self.files.main_out = "optimization.out"
            self.files.wavefunc_in = "optimization-i"
            self.files.wavefunc_out = "optimization-o"
            self.files.tmp = "tmp"
            #self.files.to_files(fout, self.dataset[0].system)
            fout.write(self.files.to_string(system=self.dataset[0].system))
            fout.write("EOF\n")

            a = self.dataset[0].system.xyz.cell[0][0]

            fout.write("v11=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.dataset[0].system.xyz.cell[2][2])

            fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
            fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")
            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${a}-${c}\n")
                    fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}-${c}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}-${c}/optimization.in\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  cat >> relax-${a}-${c}/optimization.in<<EOF\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}-${c}/optimization.in\n")
                    fout.write("  cd relax-${a}-${c}/\n")
                    fout.write("  yhrun $PMF_ABINIT < optimization.files\n")
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  mkdir relax-${a}\n")
                    fout.write("  cp optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}/optimization.in\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}/optimization.in\n")
                    fout.write("  cd relax-${a}/\n")
                    fout.write("  yhrun $PMF_ABINIT < optimization.files\n")
                    fout.write("  cd ../\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${c}\n")
                    fout.write("  cp optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${c}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${c}/optimization.in\n")
                    fout.write("  cat >> relax-${c}/optimization.in<<EOF\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("${v21} ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${c}/optimization.in\n")
                    fout.write("  cd relax-${c}/\n")
                    fout.write("  yhrun $PMF_ABINIT < optimization.files\n")
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass

        # gen local bash script
        with open("opt-hexagonal.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")

            fout.write("cat > optimization.in<<EOF\n")
            #self.dataset[0].to_input(fout)
            fout.write(self.dataset[0].to_string())
            fout.write("EOF\n")
            fout.write("cat > optimization.files<<EOF\n")
            #self.files.name = "optimization.files"
            self.files.main_in = "optimization.in"
            self.files.main_out = "optimization.out"
            self.files.wavefunc_in = "optimization-i"
            self.files.wavefunc_out = "optimization-o"
            self.files.tmp = "tmp"
            #self.files.to_files(fout, self.dataset[0].system)
            fout.write(self.files.to_string(system=self.dataset[0].system))
            fout.write("EOF\n")

            a = self.dataset[0].system.xyz.cell[0][0]

            fout.write("v11=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.dataset[0].system.xyz.cell[2][2])

            fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
            fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")
            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${a}-${c}\n")
                    fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}-${c}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}-${c}/optimization.in\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  cat >> relax-${a}-${c}/optimization.in<<EOF\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}-${c}/optimization.in\n")
                    fout.write("  cd relax-${a}-${c}/\n")
                    fout.write("  %s $PMF_ABINIT < optimization.files\n" % self.run_params["mpi"])
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  mkdir relax-${a}\n")
                    fout.write("  cp optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}/optimization.in\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}/optimization.in\n")
                    fout.write("  cd relax-${a}/\n")
                    fout.write("  %s $PMF_ABINIT < optimization.files\n" % self.run_params["mpi"])
                    fout.write("  cd ../\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${c}\n")
                    fout.write("  cp optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${c}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${c}/optimization.in\n")
                    fout.write("  cat >> relax-${c}/optimization.in<<EOF\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("${v21} ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${c}/optimization.in\n")
                    fout.write("  cd relax-${c}/\n")
                    fout.write("  %s $PMF_ABINIT < optimization.files\n" % self.run_params["mpi"])
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
                    fout.write("  energy=`cat ../relax-${a}-%{c}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 1`\n" % self.files.main_out)
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${c} ${energy}\n")
                    fout.write("EOF\n")
                    fout.write("done\n")
                    fout.write("done\n")

                    fout.write("cat > energy-latconst.gp<<EOF\n")
                    fout.write("set term gif\n")
                    fout.write("set output 'energy-latconst.gif'\n")
                    fout.write("set title Energy Latconst\n")
                    fout.write("set xlabel 'latconst(a)'\n")
                    fout.write("set ylabel 'latconst(c)'\n")
                    fout.write("set zlabel 'Energy'\n")
                    fout.write("splot 'energy-latconst.data'\n")
                    fout.write("EOF\n")
                else:
                    fout.write("  energy=`cat ../relax-${a}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 1`\n" % self.files.main_out)
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${energy}\n")
                    fout.write("EOF\n")
                    fout.write("done\n")

                    fout.write("cat > energy-latconst.gp<<EOF\n")
                    fout.write("set term gif\n")
                    fout.write("set output 'energy-latconst.gif'\n")
                    fout.write("set title Energy Latconst\n")
                    fout.write("set xlabel 'latconst(a)'\n")
                    fout.write("set ylabel 'Energy'\n")
                    fout.write("plot 'energy-latconst.data' w l\n")
                    fout.write("EOF\n")
            else:
                if nc >= 2:
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  energy=`cat ../relax-${c}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 1`\n" % self.files.main_out)
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${c} ${energy}\n")
                    fout.write("EOF\n")
                    fout.write("done\n")
                    fout.write("cat > energy-latconst.gp<<EOF\n")
                    fout.write("set term gif\n")
                    fout.write("set output 'energy-latconst.gif'\n")
                    fout.write("set title Energy Latconst\n")
                    fout.write("set xlabel 'latconst(c)'\n")
                    fout.write("set ylabel 'Energy'\n")
                    fout.write("plot 'energy-latconst.data' w l\n")
                    fout.write("EOF\n")
                else:
                    # nothing to do
                    pass
        #os.system("cd post-processing; bash get_energy.sh; cd ../")
        os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.system("bash %s" % "opt-hexagonal.sh")
            os.chdir("../")

        server_handle(auto=auto, directory=directory, jobfilebase="opt-hexagonal", server=self.run_params["server"])

    def tetragonal(self, directory="tmp-abinit-opt-tetragonal", runopt="genrun", auto=0, na=10, stepa=0.05, nc=10, stepc=0.05):
        """
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        os.system("cp *.psp8 %s/" % directory)
        os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
        os.system("cp %s %s/" % (self.dataset[0].system.xyz.file, directory))


        #task.poscar.to_poscar(os.path.join(directory, "POSCAR"))

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
            fout.write("cat > optimization.in<<EOF\n")
            #self.dataset[0].to_input(fout)
            fout.write(self.dataset[0].to_string())
            fout.write("EOF\n")
            fout.write("cat > optimization.files<<EOF\n")
            #self.files.name = "optimization.files"
            self.files.main_in = "optimization.in"
            self.files.main_out = "optimization.out"
            self.files.wavefunc_in = "optimization-i"
            self.files.wavefunc_out = "optimization-o"
            self.files.tmp = "tmp"
            #self.files.to_files(fout, self.dataset[0].system)
            fout.write(self.files.to_string(system=self.dataset[0].system))
            fout.write("EOF\n")

            a = self.dataset[0].system.xyz.cell[0][0]

            fout.write("v11=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.dataset[0].system.xyz.cell[2][2])

            fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
            fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${a}-${c}\n")
                    fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}-${c}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}-${c}/optimization.in\n")
                    fout.write("  cat >> relax-${a}-${c}/optimization.in<<EOF\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}-${c}/optimization.in\n")
                    fout.write("  cd relax-${a}-${c}/\n")
                    fout.write("  yhrun $PMF_ABINIT < optimization.files\n")
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  mkdir relax-${a}\n")
                    fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}/optimization.in\n")
                    fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}/optimization.in\n")
                    fout.write("  cd relax-${a}/\n")
                    fout.write("  yhrun $PMF_ABINIT < optimization.files\n")
                    fout.write("  cd ../\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${c}\n")
                    fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${c}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${c}/optimization.in\n")
                    fout.write("  cat >> relax-${c}/optimization.in<<EOF\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("0.000000 ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${c}/optimization.in\n")
                    fout.write("  cd relax-${c}/\n")
                    fout.write("  yhrun $PMF_ABINIT < optimization.files\n")
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
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            fout.write("cat > optimization.in<<EOF\n")
            #self.dataset[0].to_input(fout)
            fout.write(self.dataset[0].to_string())
            fout.write("EOF\n")
            fout.write("cat > optimization.files<<EOF\n")
            #self.files.name = "optimization.files"
            self.files.main_in = "optimization.in"
            self.files.main_out = "optimization.out"
            self.files.wavefunc_in = "optimization-i"
            self.files.wavefunc_out = "optimization-o"
            self.files.tmp = "tmp"
            #self.files.to_files(fout, self.dataset[0].system)
            fout.write(self.files.to_string(system=self.dataset[0].system))
            fout.write("EOF\n")

            a = self.dataset[0].system.xyz.cell[0][0]

            fout.write("v11=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.dataset[0].system.xyz.cell[2][2])

            fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
            fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${a}-${c}\n")
                    fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}-${c}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}-${c}/optimization.in\n")
                    fout.write("  cat >> relax-${a}-${c}/optimization.in<<EOF\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}-${c}/optimization.in\n")
                    fout.write("  cd relax-${a}-${c}/\n")
                    fout.write("  yhrun $PMF_ABINIT < optimization.files\n")
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  mkdir relax-${a}\n")
                    fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}/optimization.in\n")
                    fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}/optimization.in\n")
                    fout.write("  cd relax-${a}/\n")
                    fout.write("  yhrun $PMF_ABINIT < optimization.files\n")
                    fout.write("  cd ../\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${c}\n")
                    fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${c}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${c}/optimization.in\n")
                    fout.write("  cat >> relax-${c}/optimization.in<<EOF\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("0.000000 ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${c}/optimization.in\n")
                    fout.write("  cd relax-${c}/\n")
                    fout.write("  yhrun $PMF_ABINIT < optimization.files\n")
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass

        # gen local bash script
        with open("opt-tetragonal.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("cat > optimization.in<<EOF\n")
            #self.dataset[0].to_input(fout)
            fout.write(self.dataset[0].to_string())
            fout.write("EOF\n")
            fout.write("cat > optimization.files<<EOF\n")
            #self.files.name = "optimization.files"
            self.files.main_in = "optimization.in"
            self.files.main_out = "optimization.out"
            self.files.wavefunc_in = "optimization-i"
            self.files.wavefunc_out = "optimization-o"
            self.files.tmp = "tmp"
            #self.files.to_files(fout, self.dataset[0].system)
            fout.write(self.files.to_string(system=self.dataset[0].system))
            fout.write("EOF\n")

            a = self.dataset[0].system.xyz.cell[0][0]

            fout.write("v11=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.dataset[0].system.xyz.cell[2][2])

            fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
            fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${a}-${c}\n")
                    fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}-${c}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}-${c}/optimization.in\n")
                    fout.write("  cat >> relax-${a}-${c}/optimization.in<<EOF\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}-${c}/optimization.in\n")
                    fout.write("  cd relax-${a}-${c}/\n")
                    fout.write("  %s $PMF_ABINIT < optimization.files\n" % self.run_params["mpi"])
                    fout.write("  cd ../\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  mkdir relax-${a}\n")
                    fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}/optimization.in\n")
                    fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}/optimization.in\n")
                    fout.write("  cd relax-${a}/\n")
                    fout.write("  %s $PMF_ABINIT < optimization.files\n" % self.run_params["mpi"])
                    fout.write("  cd ../\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  mkdir relax-${c}\n")
                    fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${c}/\n")
                    fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${c}/optimization.in\n")
                    fout.write("  cat >> relax-${c}/optimization.in<<EOF\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("0.000000 ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${c}/optimization.in\n")
                    fout.write("  cd relax-${c}/\n")
                    fout.write("  %s $PMF_ABINIT < optimization.files\n" % self.run_params["mpi"])
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
                    fout.write("  energy=`cat ../relax-${a}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 1`\n" % self.files.main_out)
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${c} ${energy:32:-36}\n")
                    fout.write("EOF\n")
                    fout.write("done\n")
                    fout.write("done\n")

                    fout.write("cat > energy-latconst.gp<<EOF\n")
                    fout.write("set term gif\n")
                    fout.write("set output 'energy-latconst.gif'\n")
                    fout.write("set title Energy Latconst\n")
                    fout.write("set xlabel 'latconst(a)'\n")
                    fout.write("set ylabel 'latconst(c)'\n")
                    fout.write("set zlabel 'Energy'\n")
                    fout.write("splot 'energy-latconst.data'\n")
                    fout.write("EOF\n")
                else:
                    fout.write("  energy=`cat ../relax-${a}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 1`\n" % self.files.main_out)
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${energy:32:-36}\n")
                    fout.write("EOF\n")
                    fout.write("done\n")

                    fout.write("cat > energy-latconst.gp<<EOF\n")
                    fout.write("set term gif\n")
                    fout.write("set output 'energy-latconst.gif'\n")
                    fout.write("set title Energy Latconst\n")
                    fout.write("set xlabel 'latconst(a)'\n")
                    fout.write("set ylabel 'Energy'\n")
                    fout.write("plot 'energy-latconst.data' w l\n")
                    fout.write("EOF\n")
            else:
                if nc >= 2:
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  energy=`cat ../relax-${a}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 1`\n" % self.files.main_out)
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${c} ${energy:32:-36}\n")
                    fout.write("EOF\n")
                    fout.write("done\n")
                    fout.write("cat > energy-latconst.gp<<EOF\n")
                    fout.write("set term gif\n")
                    fout.write("set output 'energy-latconst.gif'\n")
                    fout.write("set title Energy Latconst\n")
                    fout.write("set xlabel 'latconst(c)'\n")
                    fout.write("set ylabel 'Energy'\n")
                    fout.write("plot 'energy-latconst.data' w l\n")
                    fout.write("EOF\n")
                else:
                    # nothing to do
                    pass
        #os.system("cd post-processing; bash get_energy.sh; cd ../")
        os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.system("bash %s" % "opt-tetragonal.sh")
            os.chdir("../")

        server_handle(auto=auto, directory=directory, jobfilebase="opt-tetragonal", server=self.run_params["server"])
