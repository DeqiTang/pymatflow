#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import numpy as np

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
            self.gen_pbs(directory=directory, script="optimization.pbs", cmd="$PMF_ABINIT", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
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
            self.gen_pbs(directory=directory, script="optimization.pbs", cmd="$PMF_ABINIT", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
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
            self.dataset[0].system.coordtype = "reduced"
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


            a = np.sqrt(self.dataset[0].system.xyz.cell[0][0]**2+self.dataset[0].system.xyz.cell[0][1]**2+self.dataset[0].system.xyz.cell[0][2]**2)
            b = np.sqrt(self.dataset[0].system.xyz.cell[1][0]**2+self.dataset[0].system.xyz.cell[1][1]**2+self.dataset[0].system.xyz.cell[1][2]**2)
            c = np.sqrt(self.dataset[0].system.xyz.cell[2][0]**2+self.dataset[0].system.xyz.cell[2][1]**2+self.dataset[0].system.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.dataset[0].system.xyz.cell[2][2])

            fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
            fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")
            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  mkdir relax-${a}\n")
            fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}/\n")
            fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}/optimization.in\n")
            fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
            fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec31=`echo \"scale=6; result=${c1} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec32=`echo \"scale=6; result=${c2} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
            fout.write("  vec33=`echo \"scale=6; result=${c3} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                      
            fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
            fout.write("${vec11} ${vec12} ${vec13}\n")
            fout.write("${vec21} ${vec22} ${vec23}\n")
            fout.write("${vec31} ${vec32} ${vec33}\n")
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
            self.dataset[0].system.coordtype = "reduced"            
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


            a = np.sqrt(self.dataset[0].system.xyz.cell[0][0]**2+self.dataset[0].system.xyz.cell[0][1]**2+self.dataset[0].system.xyz.cell[0][2]**2)
            b = np.sqrt(self.dataset[0].system.xyz.cell[1][0]**2+self.dataset[0].system.xyz.cell[1][1]**2+self.dataset[0].system.xyz.cell[1][2]**2)
            c = np.sqrt(self.dataset[0].system.xyz.cell[2][0]**2+self.dataset[0].system.xyz.cell[2][1]**2+self.dataset[0].system.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.dataset[0].system.xyz.cell[2][2])

            fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
            fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")
            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  mkdir relax-${a}\n")
            fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}/\n")
            fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}/optimization.in\n")
            fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
            fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec31=`echo \"scale=6; result=${c1} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec32=`echo \"scale=6; result=${c2} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
            fout.write("  vec33=`echo \"scale=6; result=${c3} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                          
            fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
            fout.write("${vec11} ${vec12} ${vec13}\n")
            fout.write("${vec21} ${vec22} ${vec23}\n")
            fout.write("${vec31} ${vec32} ${vec33}\n")
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
            self.dataset[0].system.coordtype = "reduced"            
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


            a = np.sqrt(self.dataset[0].system.xyz.cell[0][0]**2+self.dataset[0].system.xyz.cell[0][1]**2+self.dataset[0].system.xyz.cell[0][2]**2)
            b = np.sqrt(self.dataset[0].system.xyz.cell[1][0]**2+self.dataset[0].system.xyz.cell[1][1]**2+self.dataset[0].system.xyz.cell[1][2]**2)
            c = np.sqrt(self.dataset[0].system.xyz.cell[2][0]**2+self.dataset[0].system.xyz.cell[2][1]**2+self.dataset[0].system.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.dataset[0].system.xyz.cell[2][2])

            fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
            fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")
            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  mkdir relax-${a}\n")
            fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}/\n")
            fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}/optimization.in\n")
            fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
            fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec31=`echo \"scale=6; result=${c1} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec32=`echo \"scale=6; result=${c2} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
            fout.write("  vec33=`echo \"scale=6; result=${c3} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                          
            fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
            fout.write("${vec11} ${vec12} ${vec13}\n")
            fout.write("${vec21} ${vec22} ${vec23}\n")
            fout.write("${vec31} ${vec32} ${vec33}\n")
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
            fout.write("  energy=`cat ../relax-${a}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 2`\n" % self.files.main_out)
            fout.write("  cat >> energy-latconst.data <<EOF\n")
            fout.write("${a} ${energy}\n")
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
            fout.write("gnuplot energy-latconst.gp\n")

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
            self.dataset[0].system.coordtype = "reduced"            
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

            a = np.sqrt(self.dataset[0].system.xyz.cell[0][0]**2+self.dataset[0].system.xyz.cell[0][1]**2+self.dataset[0].system.xyz.cell[0][2]**2)
            b = np.sqrt(self.dataset[0].system.xyz.cell[1][0]**2+self.dataset[0].system.xyz.cell[1][1]**2+self.dataset[0].system.xyz.cell[1][2]**2)
            c = np.sqrt(self.dataset[0].system.xyz.cell[2][0]**2+self.dataset[0].system.xyz.cell[2][1]**2+self.dataset[0].system.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.dataset[0].system.xyz.cell[2][2])

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
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                  
                    fout.write("  cat >> relax-${a}-${c}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
                    fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                            
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                            
                    fout.write("  cat >> relax-${c}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
            self.dataset[0].system.coordtype = "reduced"            
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

            a = np.sqrt(self.dataset[0].system.xyz.cell[0][0]**2+self.dataset[0].system.xyz.cell[0][1]**2+self.dataset[0].system.xyz.cell[0][2]**2)
            b = np.sqrt(self.dataset[0].system.xyz.cell[1][0]**2+self.dataset[0].system.xyz.cell[1][1]**2+self.dataset[0].system.xyz.cell[1][2]**2)
            c = np.sqrt(self.dataset[0].system.xyz.cell[2][0]**2+self.dataset[0].system.xyz.cell[2][1]**2+self.dataset[0].system.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.dataset[0].system.xyz.cell[2][2])

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
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                            
                    fout.write("  cat >> relax-${a}-${c}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                            
                    fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")                    
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")        
                    fout.write("  cat >> relax-${c}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
            self.dataset[0].system.coordtype = "reduced"            
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

            a = np.sqrt(self.dataset[0].system.xyz.cell[0][0]**2+self.dataset[0].system.xyz.cell[0][1]**2+self.dataset[0].system.xyz.cell[0][2]**2)
            b = np.sqrt(self.dataset[0].system.xyz.cell[1][0]**2+self.dataset[0].system.xyz.cell[1][1]**2+self.dataset[0].system.xyz.cell[1][2]**2)
            c = np.sqrt(self.dataset[0].system.xyz.cell[2][0]**2+self.dataset[0].system.xyz.cell[2][1]**2+self.dataset[0].system.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.dataset[0].system.xyz.cell[2][2])

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
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                            
                    fout.write("  cat >> relax-${a}-${c}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                            
                    fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                            
                    fout.write("  cat >> relax-${c}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
                    fout.write("  energy=`cat ../relax-${a}-%{c}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 2`\n" % self.files.main_out)
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${c} ${energy}\n")
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
                    fout.write("gnuplot energy-latconst.gp\n")
                else:
                    fout.write("  energy=`cat ../relax-${a}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 1`\n" % self.files.main_out)
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${energy}\n")
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
                    fout.write("gnuplot energy-latconst.gp\n")
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
                    fout.write("set title 'Energy Latconst'\n")
                    fout.write("set xlabel 'latconst(c)'\n")
                    fout.write("set ylabel 'Energy'\n")
                    fout.write("plot 'energy-latconst.data' w l\n")
                    fout.write("EOF\n")
                    fout.write("gnuplot energy-latconst.gp\n")
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
            self.dataset[0].system.coordtype = "reduced"            
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

            a = np.sqrt(self.dataset[0].system.xyz.cell[0][0]**2+self.dataset[0].system.xyz.cell[0][1]**2+self.dataset[0].system.xyz.cell[0][2]**2)
            b = np.sqrt(self.dataset[0].system.xyz.cell[1][0]**2+self.dataset[0].system.xyz.cell[1][1]**2+self.dataset[0].system.xyz.cell[1][2]**2)
            c = np.sqrt(self.dataset[0].system.xyz.cell[2][0]**2+self.dataset[0].system.xyz.cell[2][1]**2+self.dataset[0].system.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.dataset[0].system.xyz.cell[2][2])

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
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                            
                    fout.write("  cat >> relax-${a}-${c}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")        
                    fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                            
                    fout.write("  cat >> relax-${c}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
            self.dataset[0].system.coordtype = "reduced"            
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

            a = np.sqrt(self.dataset[0].system.xyz.cell[0][0]**2+self.dataset[0].system.xyz.cell[0][1]**2+self.dataset[0].system.xyz.cell[0][2]**2)
            b = np.sqrt(self.dataset[0].system.xyz.cell[1][0]**2+self.dataset[0].system.xyz.cell[1][1]**2+self.dataset[0].system.xyz.cell[1][2]**2)
            c = np.sqrt(self.dataset[0].system.xyz.cell[2][0]**2+self.dataset[0].system.xyz.cell[2][1]**2+self.dataset[0].system.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.dataset[0].system.xyz.cell[2][2])

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
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")        
                    fout.write("  cat >> relax-${a}-${c}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")        
                    fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")        
                    fout.write("  cat >> relax-${c}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
            self.dataset[0].system.coordtype = "reduced"            
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

            a = np.sqrt(self.dataset[0].system.xyz.cell[0][0]**2+self.dataset[0].system.xyz.cell[0][1]**2+self.dataset[0].system.xyz.cell[0][2]**2)
            b = np.sqrt(self.dataset[0].system.xyz.cell[1][0]**2+self.dataset[0].system.xyz.cell[1][1]**2+self.dataset[0].system.xyz.cell[1][2]**2)
            c = np.sqrt(self.dataset[0].system.xyz.cell[2][0]**2+self.dataset[0].system.xyz.cell[2][1]**2+self.dataset[0].system.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.dataset[0].system.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.dataset[0].system.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.dataset[0].system.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.dataset[0].system.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.dataset[0].system.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.dataset[0].system.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.dataset[0].system.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.dataset[0].system.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.dataset[0].system.xyz.cell[2][2])

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
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")        
                    fout.write("  cat >> relax-${a}-${c}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")        
                    fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")        
                    fout.write("  cat >> relax-${c}/optimization.in<<EOF\n")
                    fout.write("${vec11} ${vec12} ${vec13}\n")
                    fout.write("${vec21} ${vec22} ${vec23}\n")
                    fout.write("${vec31} ${vec32} ${vec33}\n")
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
                    fout.write("  energy=`cat ../relax-${a}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 2`\n" % self.files.main_out)
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${c} ${energy:32:-36}\n")
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
                    fout.write("gnuplot energy-latconst.gp\n")
                else:
                    fout.write("  energy=`cat ../relax-${a}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 2`\n" % self.files.main_out)
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${energy:32:-36}\n")
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
                    fout.write("gnuplot energy-latconst.gp\n")
            else:
                if nc >= 2:
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  energy=`cat ../relax-${a}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 2`\n" % self.files.main_out)
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${c} ${energy:32:-36}\n")
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
                    fout.write("gnuplot energy-latconst.gp\n")
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



    def abc(self, directory="tmp-qe-opt-abc", runopt="gen", auto=0, range_a=[-0.1, 0.1, 0.01], range_b=[-0.1, 0.1, 0.01], range_c=[-0.1, 0.1, 0.01]):
        """
        """

        na = len(np.arange(range_a[0], range_a[1], range_a[2]))
        nb = len(np.arange(range_b[0], range_b[1], range_b[2]))
        nc = len(np.arange(range_c[0], range_c[1], range_c[2]))
        

        if self.batch_a == None:
            # namely all in one batch
            self.batch_a = na
        else:
            pass    

        if self.batch_b == None:
            # namely all in one batch
            self.batch_b = nb
        else:
            pass
        
        if self.batch_c == None:
            # namely all in one batch
            self.batch_c = nc
        else:
            pass

        if na % self.batch_a == 0:
            n_batch_a = int(na / self.batch_a)
        else:
            n_batch_a = int(na / self.batch_a) + 1

        if nb % self.batch_b == 0:
            n_batch_b = int(nb / self.batch_b)
        else:
            n_batch_b = int(nb / self.batch_b) + 1

        if nc % self.batch_c == 0:
            n_batch_c = int(nc / self.batch_c)
        else:
            n_batch_c = int(nc / self.batch_c) + 1        
        #

        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        os.system("cp *.psp8 %s/" % directory)
        os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
        os.system("cp %s %s/" % (self.dataset[0].system.xyz.file, directory))

        
        os.chdir(directory)

        for i_batch_a in range(n_batch_a):
            for i_batch_b in range(n_batch_b):
                for i_batch_c in range(n_batch_c):
                    
                    range_a_start = range_a[0] + i_batch_a * self.batch_a * range_a[2]
                    range_a_end = range_a[0] + (i_batch_a+1) * self.batch_a * range_a[2] - range_a[2] / 2
                    # - range_a[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_a_end  > range_a[1]:
                        range_a_end = range_a[1]

                    range_b_start = range_b[0] + i_batch_b * self.batch_b * range_b[2]
                    range_b_end = range_b[0] + (i_batch_b+1) * self.batch_b * range_b[2] - range_b[2] / 2
                    # - range_b[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_b_end  > range_b[1]:
                        range_b_end = range_b[1]
                        
                    range_c_start = range_c[0] + i_batch_c * self.batch_c * range_c[2]
                    range_c_end = range_c[0] + (i_batch_c+1) * self.batch_c * range_c[2] - range_c[2] / 2
                    # - range_c[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_c_end  > range_c[1]:
                        range_c_end = range_c[1]
                        
                    a = np.sqrt(self.dataset[0].system.xyz.cell[0][0]**2+self.dataset[0].system.xyz.cell[0][1]**2+self.dataset[0].system.xyz.cell[0][2]**2)
                    b = np.sqrt(self.dataset[0].system.xyz.cell[1][0]**2+self.dataset[0].system.xyz.cell[1][1]**2+self.dataset[0].system.xyz.cell[1][2]**2)
                    c = np.sqrt(self.dataset[0].system.xyz.cell[2][0]**2+self.dataset[0].system.xyz.cell[2][1]**2+self.dataset[0].system.xyz.cell[2][2]**2)
                        
                            
                    # gen llhpc script
                    with open("opt-abc-%d-%d-%d.slurm" % (i_batch_a, i_batch_b, i_batch_c), 'w') as fout:
                        fout.write("#!/bin/bash\n")
                        fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                        fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                        fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                        fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                        fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                        fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                        fout.write("cat > optimization.in<<EOF\n")
                        #self.dataset[0].to_input(fout)
                        self.dataset[0].system.coordtype = "reduced"            
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


                        fout.write("a_in=%f\n" % a)
                        fout.write("b_in=%f\n" % b)
                        fout.write("c_in=%f\n" % c)

                        fout.write("a1=%f\n" % self.dataset[0].system.xyz.cell[0][0])
                        fout.write("a2=%f\n" % self.dataset[0].system.xyz.cell[0][1])
                        fout.write("a3=%f\n" % self.dataset[0].system.xyz.cell[0][2])
                        fout.write("b1=%f\n" % self.dataset[0].system.xyz.cell[1][0])
                        fout.write("b2=%f\n" % self.dataset[0].system.xyz.cell[1][1])
                        fout.write("b3=%f\n" % self.dataset[0].system.xyz.cell[1][2])
                        fout.write("c1=%f\n" % self.dataset[0].system.xyz.cell[2][0])
                        fout.write("c2=%f\n" % self.dataset[0].system.xyz.cell[2][1])
                        fout.write("c3=%f\n" % self.dataset[0].system.xyz.cell[2][2])

                        fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
                        fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")


                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        fout.write("for b in `seq -w %f %f %f`\n" % (b+range_b_start, range_b[2], b+range_b_end))
                        fout.write("do\n")
                        fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                        fout.write("do\n")
                        fout.write("  mkdir relax-${a}-${b}-${c}\n")
                        fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}-${b}-${c}/\n")
                        fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}-${b}-${c}/optimization.in\n")
                        fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                   
                        fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec21=`echo \"scale=6; result=${b1} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec22=`echo \"scale=6; result=${b2} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec23=`echo \"scale=6; result=${b3} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                       
                        fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                         
                        fout.write("  cat >> relax-${a}-${b}-${c}/optimization.in<<EOF\n")
                        fout.write("${vec11} ${vec12} ${vec13}\n")
                        fout.write("${vec21} ${vec22} ${vec23}\n")
                        fout.write("${vec31} ${vec32} ${vec33}\n")
                        fout.write("EOF\n")
                        fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}-${b}-${c}/optimization.in\n")
                        fout.write("  cd relax-${a}-${b}-${c}/\n")
                        fout.write("  yhrun $PMF_ABINIT < optimization.files\n")
                        fout.write("  cd ../\n")
                        fout.write("done\n")
                        fout.write("done\n")
                        fout.write("done\n")
                    

                    

                    # gen pbs script
                    with open("opt-abc-%d-%d-%d.pbs" % (i_batch_a, i_batch_b, i_batch_c), 'w') as fout:
                        fout.write("#!/bin/bash\n")
                        fout.write("#PBS -N %s-%d-%d-%d\n" % (self.run_params["jobname"], i_batch_a, i_batch_b, i_batch_c))
                        fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                        if "queue" in self.run_params and self.run_params["queue"] != None:
                            fout.write("#PBS -q %s\n" %self.run_params["queue"])            
                        fout.write("\n")
                        fout.write("cd $PBS_O_WORKDIR\n")
                        fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")

                        fout.write("cat > optimization.in<<EOF\n")
                        #self.dataset[0].to_input(fout)
                        self.dataset[0].system.coordtype = "reduced"            
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


                        fout.write("a_in=%f\n" % a)
                        fout.write("b_in=%f\n" % b)
                        fout.write("c_in=%f\n" % c)

                        fout.write("a1=%f\n" % self.dataset[0].system.xyz.cell[0][0])
                        fout.write("a2=%f\n" % self.dataset[0].system.xyz.cell[0][1])
                        fout.write("a3=%f\n" % self.dataset[0].system.xyz.cell[0][2])
                        fout.write("b1=%f\n" % self.dataset[0].system.xyz.cell[1][0])
                        fout.write("b2=%f\n" % self.dataset[0].system.xyz.cell[1][1])
                        fout.write("b3=%f\n" % self.dataset[0].system.xyz.cell[1][2])
                        fout.write("c1=%f\n" % self.dataset[0].system.xyz.cell[2][0])
                        fout.write("c2=%f\n" % self.dataset[0].system.xyz.cell[2][1])
                        fout.write("c3=%f\n" % self.dataset[0].system.xyz.cell[2][2])

                        fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
                        fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")


                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        fout.write("for b in `seq -w %f %f %f`\n" % (b+range_b_start, range_b[2], b+range_b_end))
                        fout.write("do\n")
                        fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                        fout.write("do\n")
                        fout.write("  mkdir relax-${a}-${b}-${c}\n")
                        fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}-${b}-${c}/\n")
                        fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}-${b}-${c}/optimization.in\n")
                        fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                   
                        fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec21=`echo \"scale=6; result=${b1} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec22=`echo \"scale=6; result=${b2} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec23=`echo \"scale=6; result=${b3} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                       
                        fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                         
                        fout.write("  cat >> relax-${a}-${b}-${c}/optimization.in<<EOF\n")
                        fout.write("${vec11} ${vec12} ${vec13}\n")
                        fout.write("${vec21} ${vec22} ${vec23}\n")
                        fout.write("${vec31} ${vec32} ${vec33}\n")
                        fout.write("EOF\n")
                        fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}-${b}-${c}/optimization.in\n")
                        fout.write("  cd relax-${a}-${b}-${c}/\n")
                        fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_ABINIT < optimization.files\n")
                        fout.write("  cd ../\n")
                        fout.write("done\n")
                        fout.write("done\n")
                        fout.write("done\n")


                    # gen local bash script
                    with open("opt-abc-%d-%d-%d.sh" % (i_batch_a, i_batch_b, i_batch_c), 'w') as fout:
                        fout.write("#!/bin/bash\n")

                        fout.write("cat > optimization.in<<EOF\n")
                        #self.dataset[0].to_input(fout)
                        self.dataset[0].system.coordtype = "reduced"            
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


                        fout.write("a_in=%f\n" % a)
                        fout.write("b_in=%f\n" % b)
                        fout.write("c_in=%f\n" % c)

                        fout.write("a1=%f\n" % self.dataset[0].system.xyz.cell[0][0])
                        fout.write("a2=%f\n" % self.dataset[0].system.xyz.cell[0][1])
                        fout.write("a3=%f\n" % self.dataset[0].system.xyz.cell[0][2])
                        fout.write("b1=%f\n" % self.dataset[0].system.xyz.cell[1][0])
                        fout.write("b2=%f\n" % self.dataset[0].system.xyz.cell[1][1])
                        fout.write("b3=%f\n" % self.dataset[0].system.xyz.cell[1][2])
                        fout.write("c1=%f\n" % self.dataset[0].system.xyz.cell[2][0])
                        fout.write("c2=%f\n" % self.dataset[0].system.xyz.cell[2][1])
                        fout.write("c3=%f\n" % self.dataset[0].system.xyz.cell[2][2])

                        fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
                        fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")


                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        fout.write("for b in `seq -w %f %f %f`\n" % (b+range_b_start, range_b[2], b+range_b_end))
                        fout.write("do\n")
                        fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                        fout.write("do\n")
                        fout.write("  mkdir relax-${a}-${b}-${c}\n")
                        fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}-${b}-${c}/\n")
                        fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}-${b}-${c}/optimization.in\n")
                        fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                   
                        fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec21=`echo \"scale=6; result=${b1} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec22=`echo \"scale=6; result=${b2} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec23=`echo \"scale=6; result=${b3} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                       
                        fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                         
                        fout.write("  cat >> relax-${a}-${b}-${c}/optimization.in<<EOF\n")
                        fout.write("${vec11} ${vec12} ${vec13}\n")
                        fout.write("${vec21} ${vec22} ${vec23}\n")
                        fout.write("${vec31} ${vec32} ${vec33}\n")
                        fout.write("EOF\n")
                        fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}-${b}-${c}/optimization.in\n")
                        fout.write("  cd relax-${a}-${b}-${c}/\n")
                        fout.write("  %s $PMF_ABINIT < optimization.files\n" % (self.run_params["mpi"]))
                        fout.write("  cd ../\n")
                        fout.write("done\n")
                        fout.write("done\n")
                        fout.write("done\n")


                    # gen lsf_sz script
                    with open("opt-abc-%d-%d-%d.lsf_sz" % (i_batch_a, i_batch_b, i_batch_c), 'w') as fout:
                        fout.write("#!/bin/bash\n")
                        fout.write("APP_NAME=intelY_mid\n")
                        fout.write("NP=%d\n" % (self.run_params["nodes"]*self.run_params["ppn"]))
                        fout.write("NP_PER_NODE=%d\n" % self.run_params["ppn"])
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
                        fout.write("ndoelist=$(cat $CURDIR/nodelist | uniq | awk \'{print $1}\' | tr \'\n\' \',\')\n")


                        fout.write("cat > optimization.in<<EOF\n")
                        #self.dataset[0].to_input(fout)
                        self.dataset[0].system.coordtype = "reduced"            
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
                        
                        
                        fout.write("a_in=%f\n" % a)
                        fout.write("b_in=%f\n" % b)
                        fout.write("c_in=%f\n" % c)

                        fout.write("a1=%f\n" % self.dataset[0].system.xyz.cell[0][0])
                        fout.write("a2=%f\n" % self.dataset[0].system.xyz.cell[0][1])
                        fout.write("a3=%f\n" % self.dataset[0].system.xyz.cell[0][2])
                        fout.write("b1=%f\n" % self.dataset[0].system.xyz.cell[1][0])
                        fout.write("b2=%f\n" % self.dataset[0].system.xyz.cell[1][1])
                        fout.write("b3=%f\n" % self.dataset[0].system.xyz.cell[1][2])
                        fout.write("c1=%f\n" % self.dataset[0].system.xyz.cell[2][0])
                        fout.write("c2=%f\n" % self.dataset[0].system.xyz.cell[2][1])
                        fout.write("c3=%f\n" % self.dataset[0].system.xyz.cell[2][2])

                        fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
                        fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")


                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        fout.write("for b in `seq -w %f %f %f`\n" % (b+range_b_start, range_b[2], b+range_b_end))
                        fout.write("do\n")
                        fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                        fout.write("do\n")
                        fout.write("  mkdir relax-${a}-${b}-${c}\n")
                        fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}-${b}-${c}/\n")
                        fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}-${b}-${c}/optimization.in\n")
                        fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                   
                        fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec21=`echo \"scale=6; result=${b1} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec22=`echo \"scale=6; result=${b2} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec23=`echo \"scale=6; result=${b3} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                       
                        fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                         
                        fout.write("  cat >> relax-${a}-${b}-${c}/optimization.in<<EOF\n")
                        fout.write("${vec11} ${vec12} ${vec13}\n")
                        fout.write("${vec21} ${vec22} ${vec23}\n")
                        fout.write("${vec31} ${vec32} ${vec33}\n")
                        fout.write("EOF\n")
                        fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}-${b}-${c}/optimization.in\n")
                        fout.write("  cd relax-${a}-${b}-${c}/\n")
                        fout.write("  %s $PMF_ABINIT < optimization.files\n" % (self.run_params["mpi"]))
                        fout.write("  cd ../\n")
                        fout.write("done\n")
                        fout.write("done\n")
                        fout.write("done\n")


        # generate result analysis script
        os.system("mkdir -p post-processing")

        with open("post-processing/get_energy.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")
            # the comment
            fout.write("cat > energy-latconst.data <<EOF\n")
            fout.write("# format: a b c energy(eV)\n")
            fout.write("EOF\n")
            # end

            fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a[0], range_a[2], a+range_a[1]))
            fout.write("do\n")
            fout.write("for b in `seq -2 %f %f %f`\n" % (b+range_b[0], range_b[2], b+range_b[1]))
            fout.write("do\n")
            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c[0], range_c[2], c+range_c[1]))
            fout.write("do\n")
            fout.write("  energy=`cat ../relax-${a}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 2`\n" % self.files.main_out)
            fout.write("  cat >> energy-latconst.data <<EOF\n")
            fout.write("${a} ${b} ${c} ${energy:32:-36}\n")
            fout.write("EOF\n")
            
            #fout.write("cat > energy-latconst.gp<<EOF\n")
            #fout.write("set term gif\n")
            #fout.write("set output 'energy-latconst.gif'\n")
            #fout.write("set title 'Energy Latconst'\n")
            #fout.write("set xlabel 'latconst(a)'\n")
            #fout.write("set ylabel 'latconst(c)'\n")
            #fout.write("set zlabel 'Energy'\n")
            #fout.write("splot 'energy-latconst.data'\n")
            #fout.write("EOF\n")
            #fout.write("\n")
            #fout.write("gnuplot energy-latconst.gp")

        os.chdir("../")
        
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for i_batch_a in range(n_batch_a):
                for i_batch_b in range(n_batch_b):
                    for i_batch_c in range(n_batch_c):
                        os.system("bash opt-abc-%d-%d-%d.sh" % (i_batch_a, i_batch_b, i_batch_c))
            os.chdir("../")
        for i_batch_a in range(n_batch_a):
            for i_batch_b in range(n_batch_b):
                for i_batch_c in range(n_batch_c):
                    server_handle(auto=auto, directory=directory, jobfilebase="opt-abc-%d-%d-%d" % (i_batch_a, i_batch_b, i_batch_c), server=self.run_params["server"])
