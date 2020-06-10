"""
Optimization calculation, including GEO_OPT and CELL_OPT
"""
import numpy as np
import sys
import os
import shutil


from pymatflow.remote.server import server_handle
from pymatflow.cp2k.cp2k import cp2k
#from pymatflow.cp2k.base.glob import cp2k_glob
#from pymatflow.cp2k.base.force_eval import cp2k_force_eval
#from pymatflow.cp2k.base.motion import cp2k_motion

"""
"""

class opt_run(cp2k):
    """
    Note:
        opt_run is the calss as an agent for geometric optimization, including GEO_OPT
        and CELL_OPT.
    """
    def __init__(self):
        """
        TODO:
        """
        super().__init__()

        self.run_type = "GEO_OPT" # default is GEO_OPT, can also do CELL_OPT

        self.force_eval.basic_setting()


    def geo_opt(self, directory="tmp-cp2k-geo-opt", inpname="geo-opt.inp", output="geo-opt.out", runopt="gen", auto=0):
        """
        :param directory:
            where the calculation will happen
        :param inpname:
            input filename for the cp2k
        :param output:
            output filename for the cp2k
        """
        self.set_geo_opt()
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.motion.to_input(fout)

            # gen server job comit file
            self.gen_llhpc(directory=directory, inpname=inpname, output=output, cmd="$PMF_CP2K")
            # gen pbs server job comit file
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_CP2K", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_CP2K -in %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="geo-opt", server=self.run_params["server"])


    def cell_opt(self, directory="tmp-cp2k-cell-opt", inpname="cell-opt.inp", output="cell-opt.out", runopt="gen", auto=0):
        """
        :param directory:
            where the calculation will happen
        :param inpname:
            input filename for the cp2k
        :param output:
            output filename for the cp2k
        """
        self.set_cell_opt()
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.motion.to_input(fout)

            # gen server job comit file
            self.gen_llhpc(directory=directory, inpname=inpname, output=output, cmd="$PMF_CP2K")
            # gen pbs server job comit file
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_CP2K", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_CP2K -in %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="cell-opt", server=self.run_params["server"])

    def set_geo_opt(self):
        """
        Note:
            set basic parameters for GEO_OPT type running
        """
        self.run_type = "GEO_OPT"
        self.glob.params["RUN_TYPE"] = "GEO_OPT"
        self.motion.set_type("GEO_OPT")

    def set_cell_opt(self):
        """
        Note:
            set basic parameters for CELL_OPT type running

        Warning:
            if you are doing CELL_OPT run, you must also enable
            "STRESS_TENSOR" in FORCE_EVAL%STRESS_TENSOR
        """
        self.run_type = "CELL_OPT"
        self.glob.params["RUN_TYPE"] = "CELL_OPT"
        self.motion.set_type("CELL_OPT")
        if self.force_eval.params["STRESS_TENSOR"] is None:
            self.force_eval.params["STRESS_TENSOR"] = "ANALYTICAL" # NUMERICAL

    #
    def cubic(self, directory="tmp-cp2k-opt-cubic", runopt="gen", auto=0, na=10, stepa=0.05):
        """
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)

        shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

        #
        os.chdir(directory)

        with open("geo-opt.inp.template", 'w') as fout:
            self.glob.to_input(fout)
            self.force_eval.to_input(fout)
            self.motion.to_input(fout)

        # gen llhpc script
        with open("geo-opt-cubic.slurm", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
            fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
            fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
            fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
            fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
            fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
            fout.write("\n")
            fout.write("# get begin and end line number of the cell block in geo-opt.inp.template\n")
            fout.write("cell_block_begin=`cat geo-opt.inp.template | grep -n \'&CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("cell_block_end=`cat geo-opt.inp.template | grep -n \'&END CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("\n")

            a = np.sqrt(self.force_eval.subsys.xyz.cell[0][0]**2+self.force_eval.subsys.xyz.cell[0][1]**2+self.force_eval.subsys.xyz.cell[0][2]**2)
            b = np.sqrt(self.force_eval.subsys.xyz.cell[1][0]**2+self.force_eval.subsys.xyz.cell[1][1]**2+self.force_eval.subsys.xyz.cell[1][2]**2)
            c = np.sqrt(self.force_eval.subsys.xyz.cell[2][0]**2+self.force_eval.subsys.xyz.cell[2][1]**2+self.force_eval.subsys.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.force_eval.subsys.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.force_eval.subsys.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.force_eval.subsys.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.force_eval.subsys.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.force_eval.subsys.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.force_eval.subsys.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.force_eval.subsys.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.force_eval.subsys.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.force_eval.subsys.xyz.cell[2][2])

            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}.inp\n")
            fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
            fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec31=`echo \"scale=6; result=${c1} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec32=`echo \"scale=6; result=${c2} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
            fout.write("  vec33=`echo \"scale=6; result=${c3} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                       
            fout.write("  cat >> geo-opt-${a}.inp <<EOF\n")
            fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
            fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
            fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
            fout.write("\t\t\tPERIODIC xyz\n")
            fout.write("EOF\n")
            fout.write("  cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}.inp\n")
            fout.write("  yhrun $PMF_CP2K -inp geo-opt-${a}.inp > geo-opt-${a}.out\n")
            fout.write("done\n")

        # gen pbs script
        with open("geo-opt-cubic.pbs", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % self.run_params["jobname"])
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
            if "queue" in self.run_params and self.run_params["queue"] != None:
                fout.write("#PBS -q %s\n" %self.run_params["queue"])            
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            fout.write("\n")
            fout.write("# get begin and end line number of the cell block in geo-opt.inp.template\n")
            fout.write("cell_block_begin=`cat geo-opt.inp.template | grep -n \'&CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("cell_block_end=`cat geo-opt.inp.template | grep -n \'&END CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("\n")

            a = np.sqrt(self.force_eval.subsys.xyz.cell[0][0]**2+self.force_eval.subsys.xyz.cell[0][1]**2+self.force_eval.subsys.xyz.cell[0][2]**2)
            b = np.sqrt(self.force_eval.subsys.xyz.cell[1][0]**2+self.force_eval.subsys.xyz.cell[1][1]**2+self.force_eval.subsys.xyz.cell[1][2]**2)
            c = np.sqrt(self.force_eval.subsys.xyz.cell[2][0]**2+self.force_eval.subsys.xyz.cell[2][1]**2+self.force_eval.subsys.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.force_eval.subsys.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.force_eval.subsys.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.force_eval.subsys.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.force_eval.subsys.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.force_eval.subsys.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.force_eval.subsys.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.force_eval.subsys.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.force_eval.subsys.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.force_eval.subsys.xyz.cell[2][2])

            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}.inp\n")
            fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
            fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec31=`echo \"scale=6; result=${c1} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec32=`echo \"scale=6; result=${c2} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
            fout.write("  vec33=`echo \"scale=6; result=${c3} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                        
            fout.write("  cat >> geo-opt-${a}.inp <<EOF\n")
            fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
            fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
            fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
            fout.write("\t\t\tPERIODIC xyz\n")
            fout.write("EOF\n")
            fout.write("  cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}.inp\n")
            fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_CP2K -inp geo-opt-${a}.inp > geo-opt-${a}.out\n")
            fout.write("done\n")

        # gen local bash script
        with open("geo-opt-cubic.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")

            fout.write("\n")
            fout.write("\n")
            fout.write("# get begin and end line number of the cell block in geo-opt.inp.template\n")
            fout.write("cell_block_begin=`cat geo-opt.inp.template | grep -n \'&CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("cell_block_end=`cat geo-opt.inp.template | grep -n \'&END CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("\n")

            a = np.sqrt(self.force_eval.subsys.xyz.cell[0][0]**2+self.force_eval.subsys.xyz.cell[0][1]**2+self.force_eval.subsys.xyz.cell[0][2]**2)
            b = np.sqrt(self.force_eval.subsys.xyz.cell[1][0]**2+self.force_eval.subsys.xyz.cell[1][1]**2+self.force_eval.subsys.xyz.cell[1][2]**2)
            c = np.sqrt(self.force_eval.subsys.xyz.cell[2][0]**2+self.force_eval.subsys.xyz.cell[2][1]**2+self.force_eval.subsys.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.force_eval.subsys.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.force_eval.subsys.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.force_eval.subsys.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.force_eval.subsys.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.force_eval.subsys.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.force_eval.subsys.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.force_eval.subsys.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.force_eval.subsys.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.force_eval.subsys.xyz.cell[2][2])

            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}.inp\n")
            fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
            fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec31=`echo \"scale=6; result=${c1} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
            fout.write("  vec32=`echo \"scale=6; result=${c2} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
            fout.write("  vec33=`echo \"scale=6; result=${c3} * ${a} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                        
            fout.write("  cat >> geo-opt-${a}.inp <<EOF\n")
            fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
            fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
            fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
            fout.write("\t\t\tPERIODIC xyz\n")
            fout.write("EOF\n")
            fout.write("  cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}.inp\n")
            fout.write("  %s $PMF_CP2K -inp geo-opt-${a}.inp | tee geo-opt-${a}.out\n" % self.run_params["mpi"])
            fout.write("done\n")



        # generate result analysis script
        os.system("mkdir -p post-processing")

        with open("post-processing/get_energy.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("cat > energy-latconst.data <<EOF\n")
            fout.write("# format: a energy(Ry)\n")
            fout.write("EOF\n")
            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  energy=`cat ../geo-opt-${a}.out | grep 'ENERGY| Total FORCE_EVAL' | tail -n -1`\n")
            fout.write("  cat >> energy-latconst.data <<EOF\n")
            fout.write("${a} ${energy:48:-1}\n")
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
            os.system("bash geo-opt-cubic.sh")
            os.chdir("../")

        server_handle(auto=auto, directory=directory, jobfilebase="geo-opt-cubic", server=self.run_params["server"])

    def hexagonal(self, directory="tmp-cp2k-opt-hexagonal", runopt="gen", auto=0, na=10, nc=10, stepa=0.05, stepc=0.05):
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)

        shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

        #
        os.chdir(directory)

        with open("geo-opt.inp.template", 'w') as fout:
            self.glob.to_input(fout)
            self.force_eval.to_input(fout)
            self.motion.to_input(fout)


        # gen llhpc script
        with open("geo-opt-hexagonal.slurm", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
            fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
            fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
            fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
            fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
            fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
            fout.write("# get begin and end line number of the cell block in geo-opt.inp.template\n")
            fout.write("cell_block_begin=`cat geo-opt.inp.template | grep -n \'&CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("cell_block_end=`cat geo-opt.inp.template | grep -n \'&END CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("\n")

            a = np.sqrt(self.force_eval.subsys.xyz.cell[0][0]**2+self.force_eval.subsys.xyz.cell[0][1]**2+self.force_eval.subsys.xyz.cell[0][2]**2)
            b = np.sqrt(self.force_eval.subsys.xyz.cell[1][0]**2+self.force_eval.subsys.xyz.cell[1][1]**2+self.force_eval.subsys.xyz.cell[1][2]**2)
            c = np.sqrt(self.force_eval.subsys.xyz.cell[2][0]**2+self.force_eval.subsys.xyz.cell[2][1]**2+self.force_eval.subsys.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.force_eval.subsys.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.force_eval.subsys.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.force_eval.subsys.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.force_eval.subsys.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.force_eval.subsys.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.force_eval.subsys.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.force_eval.subsys.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.force_eval.subsys.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.force_eval.subsys.xyz.cell[2][2])

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}-${c}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")            
                    fout.write("  cat >> geo-opt-${a}-${c}.inp <<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC xyz\n")
                    fout.write("EOF\n")
                    fout.write("  cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}-${c}.inp\n")
                    fout.write("  yhrun $PMF_CP2K -inp geo-opt-${a}-${c}.inp > geo-opt-${a}-${c}.out\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                
                    fout.write("  cat >> geo-opt-${a}.inp <<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC xyz\n")
                    fout.write("EOF\n")
                    fout.write("  cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}.inp\n")
                    fout.write("  yhrun $PMF_CP2K -in geo-opt-${a}.inp > geo-opt-${a}.out\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${c}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")            
                    fout.write("  cat >> geo-opt-${c}.in<<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC xyz\n")
                    fout.write("EOF\n")
                    fout.write("  cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${c}.inp\n")
                    fout.write("  yhrun $PMF_CP2K -in geo-opt-${c}.inp > geo-opt-${c}.out\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass


        # gen pbs script
        with open("geo-opt-hexagonal.pbs", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % self.run_params["jobname"])
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
            if "queue" in self.run_params and self.run_params["queue"] != None:
                fout.write("#PBS -q %s\n" %self.run_params["queue"])            
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            fout.write("\n")
            fout.write("# get begin and end line number of the cell block in geo-opt.inp.template\n")
            fout.write("cell_block_begin=`cat geo-opt.inp.template | grep -n \'&CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("cell_block_end=`cat geo-opt.inp.template | grep -n \'&END CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("\n")

            a = np.sqrt(self.force_eval.subsys.xyz.cell[0][0]**2+self.force_eval.subsys.xyz.cell[0][1]**2+self.force_eval.subsys.xyz.cell[0][2]**2)
            b = np.sqrt(self.force_eval.subsys.xyz.cell[1][0]**2+self.force_eval.subsys.xyz.cell[1][1]**2+self.force_eval.subsys.xyz.cell[1][2]**2)
            c = np.sqrt(self.force_eval.subsys.xyz.cell[2][0]**2+self.force_eval.subsys.xyz.cell[2][1]**2+self.force_eval.subsys.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.force_eval.subsys.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.force_eval.subsys.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.force_eval.subsys.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.force_eval.subsys.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.force_eval.subsys.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.force_eval.subsys.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.force_eval.subsys.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.force_eval.subsys.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.force_eval.subsys.xyz.cell[2][2])

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}-${c}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                
                    fout.write("  cat >> geo-opt-${a}-${c}.inp <<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC xyz\n")
                    fout.write("EOF\n")
                    fout.write("  cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}-${c}.inp\n")
                    fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_CP2K -inp geo-opt-${a}-${c}.inp > geo-opt-${a}-${c}.out\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                
                    fout.write("  cat >> geo-opt-${a}.inp <<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC xyz\n")
                    fout.write("EOF\n")
                    fout.write("  cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}.inp\n")
                    fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_CP2K -in geo-opt-${a}.inp > geo-opt-${a}.out\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${c}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")            
                    fout.write("  cat >> geo-opt-${c}.in<<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC xyz\n")
                    fout.write("EOF\n")
                    fout.write("  cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${c}.inp\n")
                    fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_CP2K -in geo-opt-${c}.inp > geo-opt-${c}.out\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass

        # gen local bash script
        with open("geo-opt-hexagonal.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")

            fout.write("\n")

            fout.write("\n")
            fout.write("# get begin and end line number of the cell block in geo-opt.inp.template\n")
            fout.write("cell_block_begin=`cat geo-opt.inp.template | grep -n \'&CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("cell_block_end=`cat geo-opt.inp.template | grep -n \'&END CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("\n")

            a = np.sqrt(self.force_eval.subsys.xyz.cell[0][0]**2+self.force_eval.subsys.xyz.cell[0][1]**2+self.force_eval.subsys.xyz.cell[0][2]**2)
            b = np.sqrt(self.force_eval.subsys.xyz.cell[1][0]**2+self.force_eval.subsys.xyz.cell[1][1]**2+self.force_eval.subsys.xyz.cell[1][2]**2)
            c = np.sqrt(self.force_eval.subsys.xyz.cell[2][0]**2+self.force_eval.subsys.xyz.cell[2][1]**2+self.force_eval.subsys.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.force_eval.subsys.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.force_eval.subsys.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.force_eval.subsys.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.force_eval.subsys.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.force_eval.subsys.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.force_eval.subsys.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.force_eval.subsys.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.force_eval.subsys.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.force_eval.subsys.xyz.cell[2][2])

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}-${c}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                
                    fout.write("  cat >> geo-opt-${a}-${c}.inp <<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC xyz\n")
                    fout.write("EOF\n")
                    fout.write("  cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}-${c}.inp\n")
                    fout.write("  %s $PMF_CP2K -inp geo-opt-${a}-${c}.inp | tee geo-opt-${a}-${c}.out\n" % self.run_params["mpi"])
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                
                    fout.write("  cat >> geo-opt-${a}.inp <<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC xyz\n")
                    fout.write("EOF\n")
                    fout.write("  cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}.inp\n")
                    fout.write("  %s $PMF_CP2K -in geo-opt-${a}.inp | tee geo-opt-${a}.out\n" % self.run_params["mpi"])
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${c}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")            
                    fout.write("  cat >> geo-opt-${c}.in<<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC xyz\n")
                    fout.write("EOF\n")
                    fout.write("  cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${c}.inp\n")
                    fout.write("  %s $PMF_CP2K -in geo-opt-${c}.inp | tee geo-opt-${c}.out\n" % self.run_params["mpi"])
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
                fout.write("# format: a c energy(Ry)\n")
                fout.write("EOF\n")
            if na >= 2 and nc < 2:
                fout.write("cat > energy-latconst.data <<EOF\n")
                fout.write("# format: a energy(Ry)\n")
                fout.write("EOF\n")
            if na < 2 and nc >= 2:
                fout.write("cat > energy-latconst.data <<EOF\n")
                fout.write("# format: c energy(Ry)\n")
                fout.write("EOF\n")
            # end
            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # both a and c are optimized
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  energy=`cat ../geo-opt-${a}-${c}.out | grep 'ENERGY| Total FORCE_EVAL' | tail -n -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${c} ${energy:48:-1}\n")
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
                    fout.write("gnuplot energy-latconst.gp\n")
                else:
                    fout.write("  energy=`cat ../geo-opt-${a}.out | grep 'ENERGY| Total FORCE_EVAL' | tail -n -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${energy:48:-1}\n")
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
                    fout.write("gnuplot energy-latconst.gp\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only c is optimized
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  energy=`cat ../geo-opt-${c}.out | grep 'ENERGY| Total FORCE_EVAL' | tail -n -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${c} ${energy:48:-1}\n")
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
                    fout.write("gnuplot energy-latconst.gp\n")
                else:
                    # neither a nor c is optimized
                    pass
        #os.system("cd post-processing; bash get_energy.sh; cd ../")
        os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash geo-opt-hexagonal.sh")
            os.chdir("../")

        server_handle(auto=auto, directory=directory, jobfilebase="geo-opt-hexagonal", server=self.run_params["server"])

    def tetragonal(self, directory="tmp-cp2k-opt-tetragonal", runopt="gen", auto=0, na=10, nc=10, stepa=0.05, stepc=0.05):
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)

        shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

        #
        os.chdir(directory)

        with open("geo-opt.inp.template", 'w') as fout:
            self.glob.to_input(fout)
            self.force_eval.to_input(fout)
            self.motion.to_input(fout)


        # gen llhpc script
        with open("geo-opt-tetragonal.slurm", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
            fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
            fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
            fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
            fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
            fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
            fout.write("\n")
            fout.write("# get begin and end line number of the cell block in geo-opt.inp.template\n")
            fout.write("cell_block_begin=`cat geo-opt.inp.template | grep -n \'&CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("cell_block_end=`cat geo-opt.inp.template | grep -n \'&END CELL\' | head -n 1 | cut -d \':\' -f1\n`")
            fout.write("\n")

            a = np.sqrt(self.force_eval.subsys.xyz.cell[0][0]**2+self.force_eval.subsys.xyz.cell[0][1]**2+self.force_eval.subsys.xyz.cell[0][2]**2)
            b = np.sqrt(self.force_eval.subsys.xyz.cell[1][0]**2+self.force_eval.subsys.xyz.cell[1][1]**2+self.force_eval.subsys.xyz.cell[1][2]**2)
            c = np.sqrt(self.force_eval.subsys.xyz.cell[2][0]**2+self.force_eval.subsys.xyz.cell[2][1]**2+self.force_eval.subsys.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.force_eval.subsys.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.force_eval.subsys.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.force_eval.subsys.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.force_eval.subsys.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.force_eval.subsys.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.force_eval.subsys.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.force_eval.subsys.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.force_eval.subsys.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.force_eval.subsys.xyz.cell[2][2])

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("  for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("  do\n")
                    fout.write("    cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}-${c}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                
                    fout.write("    cat >> geo-opt-${a}-${c}.inp <<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC XYZ\n")
                    fout.write("EOF\n")
                    fout.write("    cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}-${c}.inp\n")
                    fout.write("    yhrun $PMF_CP2K -in geo-opt-${a}-${c}.inp > geo-opt-${a}-${c}.out\n")
                    fout.write("  done\n")
                else:
                    # only optimize a
                    fout.write("    cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")            
                    fout.write("    cat >> geo-opt-${a}.inp <<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC XYZ\n")
                    fout.write("EOF\n")
                    fout.write("    cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}.inp\n")
                    fout.write("    yhrun $PMF_CP2K -in geo-opt-${a}.inp > geo-opt-${a}.out\n")
                fout.write("done\n")
            else:
            # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${c}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                
                    fout.write("  cat >> geo-opt-${c}.in<<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC XYZ\n")
                    fout.write("EOF\n")
                    fout.write("    cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}.inp\n")
                    fout.write("    yhrun $PMF_CP2K -in geo-opt-${a}.inp > geo-opt-${a}.out\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass


        # gen pbs script
        with open("geo-opt-tetragonal.pbs", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % self.run_params["jobname"])
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
            if "queue" in self.run_params and self.run_params["queue"] != None:
                fout.write("#PBS -q %s\n" %self.run_params["queue"])            
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            fout.write("\n")
            fout.write("# get begin and end line number of the cell block in geo-opt.inp.template\n")
            fout.write("cell_block_begin=`cat geo-opt.inp.template | grep -n \'&CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("cell_block_end=`cat geo-opt.inp.template | grep -n \'&END CELL\' | head -n 1 | cut -d \':\' -f1\n`")
            fout.write("\n")

            a = np.sqrt(self.force_eval.subsys.xyz.cell[0][0]**2+self.force_eval.subsys.xyz.cell[0][1]**2+self.force_eval.subsys.xyz.cell[0][2]**2)
            b = np.sqrt(self.force_eval.subsys.xyz.cell[1][0]**2+self.force_eval.subsys.xyz.cell[1][1]**2+self.force_eval.subsys.xyz.cell[1][2]**2)
            c = np.sqrt(self.force_eval.subsys.xyz.cell[2][0]**2+self.force_eval.subsys.xyz.cell[2][1]**2+self.force_eval.subsys.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.force_eval.subsys.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.force_eval.subsys.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.force_eval.subsys.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.force_eval.subsys.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.force_eval.subsys.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.force_eval.subsys.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.force_eval.subsys.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.force_eval.subsys.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.force_eval.subsys.xyz.cell[2][2])

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("  for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("  do\n")
                    fout.write("    cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}-${c}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")            
                    fout.write("    cat >> geo-opt-${a}-${c}.inp <<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC XYZ\n")
                    fout.write("EOF\n")
                    fout.write("    cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}-${c}.inp\n")
                    fout.write("    mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_CP2K -in geo-opt-${a}-${c}.inp > geo-opt-${a}-${c}.out\n")
                    fout.write("  done\n")
                else:
                    # only optimize a
                    fout.write("    cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")            
                    fout.write("    cat >> geo-opt-${a}.inp <<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC XYZ\n")
                    fout.write("EOF\n")
                    fout.write("    cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}.inp\n")
                    fout.write("    mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_CP2K -in geo-opt-${a}.inp > geo-opt-${a}.out\n")
                fout.write("done\n")
            else:
            # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${c}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")            
                    fout.write("  cat >> geo-opt-${c}.in<<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC XYZ\n")
                    fout.write("EOF\n")
                    fout.write("    cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}.inp\n")
                    fout.write("    mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_CP2K -in geo-opt-${a}.inp > geo-opt-${a}.out\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass

        # gen local bash script
        with open("geo-opt-tetragonal.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")

            fout.write("\n")
            fout.write("\n")
            fout.write("# get begin and end line number of the cell block in geo-opt.inp.template\n")
            fout.write("cell_block_begin=`cat geo-opt.inp.template | grep -n \'&CELL\' | head -n 1 | cut -d \':\' -f1`\n")
            fout.write("cell_block_end=`cat geo-opt.inp.template | grep -n \'&END CELL\' | head -n 1 | cut -d \':\' -f1\n`")
            fout.write("\n")

            a = np.sqrt(self.force_eval.subsys.xyz.cell[0][0]**2+self.force_eval.subsys.xyz.cell[0][1]**2+self.force_eval.subsys.xyz.cell[0][2]**2)
            b = np.sqrt(self.force_eval.subsys.xyz.cell[1][0]**2+self.force_eval.subsys.xyz.cell[1][1]**2+self.force_eval.subsys.xyz.cell[1][2]**2)
            c = np.sqrt(self.force_eval.subsys.xyz.cell[2][0]**2+self.force_eval.subsys.xyz.cell[2][1]**2+self.force_eval.subsys.xyz.cell[2][2]**2)

            fout.write("a_in=%f\n" % a)
            fout.write("b_in=%f\n" % b)
            fout.write("c_in=%f\n" % c)

            fout.write("a1=%f\n" % self.force_eval.subsys.xyz.cell[0][0])
            fout.write("a2=%f\n" % self.force_eval.subsys.xyz.cell[0][1])
            fout.write("a3=%f\n" % self.force_eval.subsys.xyz.cell[0][2])
            fout.write("b1=%f\n" % self.force_eval.subsys.xyz.cell[1][0])
            fout.write("b2=%f\n" % self.force_eval.subsys.xyz.cell[1][1])
            fout.write("b3=%f\n" % self.force_eval.subsys.xyz.cell[1][2])
            fout.write("c1=%f\n" % self.force_eval.subsys.xyz.cell[2][0])
            fout.write("c2=%f\n" % self.force_eval.subsys.xyz.cell[2][1])
            fout.write("c3=%f\n" % self.force_eval.subsys.xyz.cell[2][2])

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("  for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("  do\n")
                    fout.write("    cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}-${c}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")            
                    fout.write("    cat >> geo-opt-${a}-${c}.inp <<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC XYZ\n")
                    fout.write("EOF\n")
                    fout.write("    cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}-${c}.inp\n")
                    fout.write("    %s $PMF_CP2K -in geo-opt-${a}-${c}.inp | tee geo-opt-${a}-${c}.out\n" % self.run_params["mpi"])
                    fout.write("  done\n")
                else:
                    # only optimize a
                    fout.write("    cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")            
                    fout.write("    cat >> geo-opt-${a}.inp <<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC XYZ\n")
                    fout.write("EOF\n")
                    fout.write("    cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}.inp\n")
                    fout.write("    %s $PMF_CP2K -in geo-opt-${a}.inp | tee geo-opt-${a}.out\n" % self.run_params["mpi"])
                fout.write("done\n")
            else:
            # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${c}.inp\n")
                    fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                    
                    fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec21=`echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec23=`echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                    fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                        
                    fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")            
                    fout.write("  cat >> geo-opt-${c}.in<<EOF\n")
                    fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                    fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                    fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                    fout.write("\t\t\tPERIODIC XYZ\n")
                    fout.write("EOF\n")
                    fout.write("    cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}.inp\n")
                    fout.write("    %s $PMF_CP2K -in geo-opt-${a}.inp | tee geo-opt-${a}.out\n" % self.run_params["mpi"])
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
                fout.write("# format: a c energy(Ry)\n")
                fout.write("EOF\n")
            if na >= 2 and nc < 2:
                fout.write("cat > energy-latconst.data <<EOF\n")
                fout.write("# format: a energy(Ry)\n")
                fout.write("EOF\n")
            if na < 2 and nc >= 2:
                fout.write("cat > energy-latconst.data <<EOF\n")
                fout.write("# format: c energy(Ry)\n")
                fout.write("EOF\n")
            # end
            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # both a and c are optimized
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("  do\n")
                    fout.write("    energy=`cat ../geo-opt-${a}-${c}.out | grep 'ENERGY| Total FORCE_EVAL' | tail -n -1`\n")
                    fout.write("    cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${c} ${energy:48:-1\n")
                    fout.write("EOF\n")
                    fout.write("  done\n")
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
                    fout.write("gnuplot energy-latconst.gp\n")
                else:
                    fout.write("  energy=`cat ../geo-opt-${a}.out | grep 'ENERGY| Total FORCE_EVAL' | tail -n -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${energy:48:-1}\n")
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
                    fout.write("gnuplot energy-latconst.gp\n")
            else:
            # a is not optimized
                if nc >= 2:
                    # only c is optimized
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  energy=`cat geo-opt-${c}.out | grep 'ENERGY| Total FORCE_EVAL' | tail -n -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${c} ${energy:48:-1}\n")
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
                    fout.write("gnuplot energy-latconst.gp\n")
                else:
                    # neither a nor c is optimized
                    pass
        #os.system("cd post-processing; bash get_energy.sh; cd ../")
        os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash geo-opt-tetragonal.sh")
            os.chdir("../")

        server_handle(auto=auto, directory=directory, jobfilebase="geo-opt-tetragonal", server=self.run_params["server"])
    
    
    def abc(self, directory="tmp-cp2k-opt-abc", runopt="gen", auto=0, range_a=[-0.1, 0.1, 0.01], range_b=[-0.1, 0.1, 0.01], range_c=[-0.1, 0.1, 0.01]):
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

        shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

        #
        os.chdir(directory)

        with open("geo-opt.inp.template", 'w') as fout:
            self.glob.to_input(fout)
            self.force_eval.to_input(fout)
            self.motion.to_input(fout)

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
                        
                    a = np.sqrt(self.force_eval.subsys.xyz.cell[0][0]**2+self.force_eval.subsys.xyz.cell[0][1]**2+self.force_eval.subsys.xyz.cell[0][2]**2)
                    b = np.sqrt(self.force_eval.subsys.xyz.cell[1][0]**2+self.force_eval.subsys.xyz.cell[1][1]**2+self.force_eval.subsys.xyz.cell[1][2]**2)
                    c = np.sqrt(self.force_eval.subsys.xyz.cell[2][0]**2+self.force_eval.subsys.xyz.cell[2][1]**2+self.force_eval.subsys.xyz.cell[2][2]**2)

                            
                    # gen llhpc script
                    with open("opt-abc-%d-%d-%d.slurm" % (i_batch_a, i_batch_b, i_batch_c), 'w') as fout:
                        fout.write("#!/bin/bash\n")
                        fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                        fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                        fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                        fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                        fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                        fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                        fout.write("\n")
                        fout.write("# get begin and end line number of the cell block in geo-opt.inp.template\n")
                        fout.write("cell_block_begin=`cat geo-opt.inp.template | grep -n \'&CELL\' | head -n 1 | cut -d \':\' -f1`\n")
                        fout.write("cell_block_end=`cat geo-opt.inp.template | grep -n \'&END CELL\' | head -n 1 | cut -d \':\' -f1\n`")
                        fout.write("\n")

                        fout.write("a_in=%f\n" % a)
                        fout.write("b_in=%f\n" % b)
                        fout.write("c_in=%f\n" % c)
                        
                        fout.write("a1=%f\n" % self.force_eval.subsys.xyz.cell[0][0])
                        fout.write("a2=%f\n" % self.force_eval.subsys.xyz.cell[0][1])
                        fout.write("a3=%f\n" % self.force_eval.subsys.xyz.cell[0][2])
                        fout.write("b1=%f\n" % self.force_eval.subsys.xyz.cell[1][0])
                        fout.write("b2=%f\n" % self.force_eval.subsys.xyz.cell[1][1])
                        fout.write("b3=%f\n" % self.force_eval.subsys.xyz.cell[1][2])
                        fout.write("c1=%f\n" % self.force_eval.subsys.xyz.cell[2][0])
                        fout.write("c2=%f\n" % self.force_eval.subsys.xyz.cell[2][1])
                        fout.write("c3=%f\n" % self.force_eval.subsys.xyz.cell[2][2])


                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        fout.write("for b in `seq -w %f %f %f`\n" % (b+range_b_start, range_b[2], b+range_b_end))
                        fout.write("do\n")
                        fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                        fout.write("do\n")
                        fout.write("    cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}-${b}-${c}.inp\n")
                        fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                   
                        fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec21=`echo \"scale=6; result=${b1} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec22=`echo \"scale=6; result=${b2} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec23=`echo \"scale=6; result=${b3} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                       
                        fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("    cat >> geo-opt-${a}-${b}-${c}.inp <<EOF\n")
                        fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                        fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                        fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                        fout.write("\t\t\tPERIODIC XYZ\n")
                        fout.write("EOF\n")
                        fout.write("    cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}-${c}-${c}.inp\n")
                        fout.write("    %s $PMF_CP2K -in geo-opt-${a}-${b}-${c}.inp | tee geo-opt-${a}-${b}-${c}.out\n" % self.run_params["mpi"])
                        fout.write("done\n")
                        fout.write("done\n")
                        fout.write("done\n")
                    

                    

                    # gen pbs script
                    with open("opt-tetragonal-%d-%d-%d.pbs" % (i_batch_a, i_batch_b, i_batch_c), 'w') as fout:
                        fout.write("#!/bin/bash\n")
                        fout.write("#PBS -N %s-%d-%d-%d\n" % (self.run_params["jobname"], i_batch_a, i_batch_b, i_batch_c))
                        fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                        if "queue" in self.run_params and self.run_params["queue"] != None:
                            fout.write("#PBS -q %s\n" %self.run_params["queue"])            
                        fout.write("\n")
                        fout.write("cd $PBS_O_WORKDIR\n")
                        fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")

                        fout.write("\n")
                        fout.write("# get begin and end line number of the cell block in geo-opt.inp.template\n")
                        fout.write("cell_block_begin=`cat geo-opt.inp.template | grep -n \'&CELL\' | head -n 1 | cut -d \':\' -f1`\n")
                        fout.write("cell_block_end=`cat geo-opt.inp.template | grep -n \'&END CELL\' | head -n 1 | cut -d \':\' -f1\n`")
                        fout.write("\n")

                        fout.write("a_in=%f\n" % a)
                        fout.write("b_in=%f\n" % b)
                        fout.write("c_in=%f\n" % c)

                        fout.write("a1=%f\n" % self.force_eval.subsys.xyz.cell[0][0])
                        fout.write("a2=%f\n" % self.force_eval.subsys.xyz.cell[0][1])
                        fout.write("a3=%f\n" % self.force_eval.subsys.xyz.cell[0][2])
                        fout.write("b1=%f\n" % self.force_eval.subsys.xyz.cell[1][0])
                        fout.write("b2=%f\n" % self.force_eval.subsys.xyz.cell[1][1])
                        fout.write("b3=%f\n" % self.force_eval.subsys.xyz.cell[1][2])
                        fout.write("c1=%f\n" % self.force_eval.subsys.xyz.cell[2][0])
                        fout.write("c2=%f\n" % self.force_eval.subsys.xyz.cell[2][1])
                        fout.write("c3=%f\n" % self.force_eval.subsys.xyz.cell[2][2])


                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        fout.write("for b in `seq -w %f %f %f`\n" % (b+range_b_start, range_b[2], b+range_b_end))
                        fout.write("do\n")
                        fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                        fout.write("do\n")
                        fout.write("    cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}-${b}-${c}.inp\n")
                        fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                   
                        fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec21=`echo \"scale=6; result=${b1} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec22=`echo \"scale=6; result=${b2} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec23=`echo \"scale=6; result=${b3} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                       
                        fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("    cat >> geo-opt-${a}-${b}-${c}.inp <<EOF\n")
                        fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                        fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                        fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                        fout.write("\t\t\tPERIODIC XYZ\n")
                        fout.write("EOF\n")
                        fout.write("    cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}-${c}-${c}.inp\n")
                        fout.write("    mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_CP2K -in geo-opt-${a}-${b}-${c}.inp | tee geo-opt-${a}-${b}-${c}.out\n")
                        fout.write("done\n")
                        fout.write("done\n")
                        fout.write("done\n")
                        


                    # gen local bash script
                    with open("opt-tetragonal.sh", 'w') as fout:
                        fout.write("#!/bin/bash\n")

                        fout.write("\n")
                        fout.write("# get begin and end line number of the cell block in geo-opt.inp.template\n")
                        fout.write("cell_block_begin=`cat geo-opt.inp.template | grep -n \'&CELL\' | head -n 1 | cut -d \':\' -f1`\n")
                        fout.write("cell_block_end=`cat geo-opt.inp.template | grep -n \'&END CELL\' | head -n 1 | cut -d \':\' -f1\n`")
                        fout.write("\n")


                        fout.write("a_in=%f\n" % a)
                        fout.write("b_in=%f\n" % b)
                        fout.write("c_in=%f\n" % c)

                        fout.write("a1=%f\n" % self.force_eval.subsys.xyz.cell[0][0])
                        fout.write("a2=%f\n" % self.force_eval.subsys.xyz.cell[0][1])
                        fout.write("a3=%f\n" % self.force_eval.subsys.xyz.cell[0][2])
                        fout.write("b1=%f\n" % self.force_eval.subsys.xyz.cell[1][0])
                        fout.write("b2=%f\n" % self.force_eval.subsys.xyz.cell[1][1])
                        fout.write("b3=%f\n" % self.force_eval.subsys.xyz.cell[1][2])
                        fout.write("c1=%f\n" % self.force_eval.subsys.xyz.cell[2][0])
                        fout.write("c2=%f\n" % self.force_eval.subsys.xyz.cell[2][1])
                        fout.write("c3=%f\n" % self.force_eval.subsys.xyz.cell[2][2])


                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        fout.write("for b in `seq -w %f %f %f`\n" % (b+range_b_start, range_b[2], b+range_b_end))
                        fout.write("do\n")
                        fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                        fout.write("do\n")
                        fout.write("    cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}-${b}-${c}.inp\n")
                        fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                   
                        fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec21=`echo \"scale=6; result=${b1} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec22=`echo \"scale=6; result=${b2} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec23=`echo \"scale=6; result=${b3} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                       
                        fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("    cat >> geo-opt-${a}-${b}-${c}.inp <<EOF\n")
                        fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                        fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                        fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                        fout.write("\t\t\tPERIODIC XYZ\n")
                        fout.write("EOF\n")
                        fout.write("    cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}-${c}-${c}.inp\n")
                        fout.write("    %s $PMF_CP2K -in geo-opt-${a}-${b}-${c}.inp | tee geo-opt-${a}-${b}-${c}.out\n" % self.run_params["mpi"])
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


                        fout.write("\n")
                        fout.write("# get begin and end line number of the cell block in geo-opt.inp.template\n")
                        fout.write("cell_block_begin=`cat geo-opt.inp.template | grep -n \'&CELL\' | head -n 1 | cut -d \':\' -f1`\n")
                        fout.write("cell_block_end=`cat geo-opt.inp.template | grep -n \'&END CELL\' | head -n 1 | cut -d \':\' -f1\n`")
                        fout.write("\n")
                        
                        fout.write("a_in=%f\n" % a)
                        fout.write("b_in=%f\n" % b)
                        fout.write("c_in=%f\n" % c)

                        fout.write("a1=%f\n" % self.force_eval.subsys.xyz.cell[0][0])
                        fout.write("a2=%f\n" % self.force_eval.subsys.xyz.cell[0][1])
                        fout.write("a3=%f\n" % self.force_eval.subsys.xyz.cell[0][2])
                        fout.write("b1=%f\n" % self.force_eval.subsys.xyz.cell[1][0])
                        fout.write("b2=%f\n" % self.force_eval.subsys.xyz.cell[1][1])
                        fout.write("b3=%f\n" % self.force_eval.subsys.xyz.cell[1][2])
                        fout.write("c1=%f\n" % self.force_eval.subsys.xyz.cell[2][0])
                        fout.write("c2=%f\n" % self.force_eval.subsys.xyz.cell[2][1])
                        fout.write("c3=%f\n" % self.force_eval.subsys.xyz.cell[2][2])

                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        fout.write("for b in `seq -w %f %f %f`\n" % (b+range_b_start, range_b[2], b+range_b_end))
                        fout.write("do\n")
                        fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                        fout.write("do\n")
                        fout.write("    cat geo-opt.inp.template | head -n +${cell_block_begin} > geo-opt-${a}-${b}-${c}.inp\n")
                        fout.write("  vec11=`echo \"scale=6; result=${a1} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec12=`echo \"scale=6; result=${a2} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                   
                        fout.write("  vec13=`echo \"scale=6; result=${a3} * ${a} / ${a_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec21=`echo \"scale=6; result=${b1} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec22=`echo \"scale=6; result=${b2} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec23=`echo \"scale=6; result=${b3} * ${b} / ${b_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec31=`echo \"scale=6; result=${c1} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("  vec32=`echo \"scale=6; result=${c2} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")                                       
                        fout.write("  vec33=`echo \"scale=6; result=${c3} * ${c} / ${c_in}; if (length(result)==scale(result) && result!=0) print 0; print result\" | bc`\n")
                        fout.write("    cat >> geo-opt-${a}-${b}-${c}.inp <<EOF\n")
                        fout.write("\t\t\tA ${vec11} ${vec12} ${vec13}\n")
                        fout.write("\t\t\tB ${vec21} ${vec22} ${vec23}\n")
                        fout.write("\t\t\tC ${vec31} ${vec32} ${vec33}\n")
                        fout.write("\t\t\tPERIODIC XYZ\n")
                        fout.write("EOF\n")
                        fout.write("    cat geo-opt.inp.template | tail -n +${cell_block_end} >> geo-opt-${a}-${c}-${c}.inp\n")
                        fout.write("    %s $PMF_CP2K -in geo-opt-${a}-${b}-${c}.inp | tee geo-opt-${a}-${b}-${c}.out\n" % self.run_params["mpi"])
                        fout.write("done\n")
                        fout.write("done\n")
                        fout.write("done\n")


        # generate result analysis script
        os.system("mkdir -p post-processing")

        with open("post-processing/get_energy.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")
            # the comment
            fout.write("cat > energy-latconst.data <<EOF\n")
            fout.write("# format: a b c energy(Ry)\n")
            fout.write("EOF\n")
            # end

            fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a[0], range_a[2], a+range_a[1]))
            fout.write("do\n")
            fout.write("for b in `seq -2 %f %f %f`\n" % (b+range_b[0], range_b[2], b+range_b[1]))
            fout.write("do\n")
            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c[0], range_c[2], c+range_c[1]))
            fout.write("do\n")
            fout.write("  energy=`cat ../geo-opt-${a}-${b}-${c}.out | grep 'ENERGY| Total FORCE_EVAL' | tail -n -1`\n")                    
            fout.write("  cat >> energy-latconst.data <<EOF\n")
            fout.write("${a} ${b} ${c} ${energy:48:-1}\n")
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
