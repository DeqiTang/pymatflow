"""
Geometric Optimization calc
"""
import os
import re
import shutil
import numpy as np
import pymatflow.base as base
from pymatflow.remote.server import server_handle
from pymatflow.qe.pwscf import PwScf

class OptRun(PwScf):
    """
    structural optimization uses both energies and forces to locate the minima
    along serach directions. usually insufficient scf convergence will lead to
    bad convergence of BFGS algorithm or even to errors. so when doing geometric
    optimization, we better set parameters to get a good scf convergece.

    when you structure is small, use a large kpoint set, or the optimization
    will not be reliable. if you structure is big enough, a small kpoint set
    will usually suffice the requirement.
    """
    def __init__(self):
        super().__init__()

        self.arts.ifstatic = False

    def relax(self, directory="tmp-qe-relax", inpname="relax.in", output="relax.out", runopt="gen", auto=0):
        """
        :param directory: a place for all the generated files
        """
        #self.set_relax()
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            #os.system("cp *.UPF %s/" % directory)
            #os.system("cp %s %s/" % (self.arts.xyz.file, directory))

            # do not copy too many files at the same time or it will be slow
            # so we do not copy all UPF files in the directory but just copy
            # those used in the calculation.
            shutil.copyfile(self.arts.xyz.file, os.path.join(directory, os.path.basename(self.arts.xyz.file)))
            #all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
            all_file = os.listdir()
            for element in self.arts.xyz.specie_labels:
                for item in all_file:
                    #if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                    #if item.split(".")[0].lower() == element.lower() or item.split("_")[0].lower() == element.lower():
                    if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE) or re.match("(%s)(_*)(upf)" % element, item, re.IGNORECASE):    
                        shutil.copyfile(item, os.path.join(directory, item))
                        break
            self.arts.pseudo.dir = os.path.abspath(directory)
            self.control.set_params({"pseudo_dir": os.path.abspath(directory)})            
            #

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.ions.to_in(fout)
                self.arts.to_in(fout)
            # gen yhbatch script
            self.gen_llhpc(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cdcloud script
            self.gen_cdcloud(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_PWX < %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="relax", server=self.run_params["server"])

    def vc_relax(self, directory="tmp-qe-vc-relax", inpname="vc-relax.in", output="vc-relax.out", runopt="gen", auto=0):
        """
        :param directory: a place for all the generated files
        """
        #self.set_vc_relax()
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            #os.system("cp *.UPF %s/" % directory)
            #os.system("cp %s %s/" % (self.arts.xyz.file, directory))

            # do not copy too many files at the same time or it will be slow
            # so we do not copy all UPF files in the directory but just copy
            # those used in the calculation.
            shutil.copyfile(self.arts.xyz.file, os.path.join(directory, os.path.basename(self.arts.xyz.file)))
            #all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
            all_file = os.listdir()
            for element in self.arts.xyz.specie_labels:
                for item in all_file:
                    #if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                    #if item.split(".")[0].lower() == element.lower() or item.split("_")[0].lower() == element.lower():
                    if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE) or re.match("(%s)(_*)(upf)" % element, item, re.IGNORECASE):    
                        shutil.copyfile(item, os.path.join(directory, item))
                        break
            #
            self.arts.pseudo.dir = os.path.abspath(directory)
            self.control.set_params({"pseudo_dir": os.path.abspath(directory)})    


            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.ions.to_in(fout)
                self.cell.to_in(fout)
                self.arts.to_in(fout)
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, cmd="$PMF_PWX", output=output, jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cdcloud script
            self.gen_cdcloud(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_PWX < %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="vc-relax", server=self.run_params["server"])

    def set_relax(self):
        self.control.calculation("relax")
        self.control.basic_setting("relax")

        self.system.basic_setting(self.arts)
        self.electrons.basic_setting()
        self.ions.basic_setting()

    def set_vc_relax(self):
        self.control.calculation("vc-relax")
        self.control.basic_setting("vc-relax")

        self.system.basic_setting(self.arts)
        self.electrons.basic_setting()
        self.ions.basic_setting()

    def cubic(self, directory="tmp-qe-relax-cubic", runopt="gen", auto=0, range_a=[-0.1, 0.101, 0.01]):
        """
        """
        
        na = len(np.arange(range_a[0], range_a[1], range_a[2]))

        if self.batch_a == None:
            # namely all in one batch
            self.batch_a = na
        else:
            pass
        
        if na % self.batch_a == 0:
            n_batch_a = int(na / self.batch_a)
        else:
            n_batch_a = int(na / self.batch_a) + 1
        #
        #         
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)

        shutil.copyfile(self.arts.xyz.file, os.path.join(directory, os.path.basename(self.arts.xyz.file)))
        #all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
        all_file = os.listdir()
        for element in self.arts.xyz.specie_labels:
            for item in all_file:
                #if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                #if item.split(".")[0].lower() == element.lower() or item.split("_")[0].lower() == element.lower():
                if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE) or re.match("(%s)(_*)(upf)" % element, item, re.IGNORECASE):    
                    shutil.copyfile(item, os.path.join(directory, item))
                    break
        self.arts.pseudo.dir = os.path.abspath(directory)
        self.control.set_params({"pseudo_dir": os.path.abspath(directory)})
        #
        os.chdir(directory)

        with open("relax.in.template", 'w') as fout:
            self.control.to_in(fout)
            self.system.to_in(fout)
            self.electrons.to_in(fout)
            self.ions.to_in(fout)

            coordtype = "crystal" # use crystal here so we could only change cell when opt cell
            fout.write("ATOMIC_SPECIES\n")
            all_file = os.listdir(self.arts.pseudo.dir)
            for element in self.arts.xyz.specie_labels:
                for item in all_file:
                    if re.match("(%s)(.*)(upf)" % (element), item, re.IGNORECASE):
                        fout.write("%s %f %s\n" % (element, base.element[element].mass, item))
                        break
            fout.write("\n")
            if coordtype == "angstrom":
                fout.write("ATOMIC_POSITIONS angstrom\n")
                if self.arts.ifstatic == True:
                    for atom in self.arts.xyz.atoms:
                        fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
                elif self.arts.ifstatic == False:
                    for atom in self.arts.xyz.atoms:
                        fout.write("%s\t%.9f\t%.9f\t%.9f" % (atom.name, atom.x, atom.y, atom.z))
                        for fix in atom.fix:
                            if fix == True:
                                fout.write("\t0")
                            elif fix == False:
                                fout.write("\t1")
                        fout.write("\n")
                else:
                    print("===============================================\n")
                    print("warning: qe.base.arts.to_in():\n")
                    print("arts.ifstatic could only be True or False\n")
                    sys.exit(1)
                fout.write("\n")
            elif coordtype == "crystal":
                # crystal namely fractional coordinate can be convert from cartesian coordinates
                # the conversion process is like transformation of presentation in quantum mechanics
                # the convmat is bulid to do the conversion
                #latcell = np.array(self.xyz.cell)
                #latcell = latcell.reshape(3, 3)
                latcell = np.array(self.arts.xyz.cell)
                convmat = np.linalg.inv(latcell.T)
                crystal_coord = np.zeros([self.arts.xyz.natom, 3])
                for i in range(self.arts.xyz.natom):
                    crystal_coord[i] = convmat.dot(np.array([self.arts.xyz.atoms[i].x, self.arts.xyz.atoms[i].y, self.arts.xyz.atoms[i].z]))
                #
                fout.write("ATOMIC_POSITIONS crystal\n")
                if self.arts.ifstatic == True:
                    for k in range(self.arts.xyz.natom):
                        fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (self.arts.xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))
                elif self.arts.ifstatic == False:
                    for k in range(self.arts.xyz.natom):
                        fout.write("%s\t%.9f\t%.9f\t%.9f" % (self.arts.xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))
                        for fix in self.arts.xyz.atoms[k].fix:
                            if fix == True:
                                fout.write("\t0")
                            elif fix == False:
                                fout.write("\t1")
                        fout.write("\n")
                else:
                    print("===============================================\n")
                    print("warning: qe.base.arts.to_in():\n")
                    print("arts.ifstatic could only be True or False\n")
                    sys.exit(1)
                fout.write("\n")
            # end crystal type ATOMIC_POSITIONS

            # writing KPOINTS to the fout
            self.arts.write_kpoints(fout)
            # =========================
            #
            # writing forces act on atoms
            if self.arts.atomic_forces_status == True:
                self.arts.write_atomic_forces(fout)
            # =========================
        for i_batch_a in range(n_batch_a):
            # gen llhpc script
            with open("relax-cubic-%d.slurm" % (i_batch_a), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s-%d\n" % (self.run_params["jobname"], i_batch_a))
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

                a = np.sqrt(self.arts.xyz.cell[0][0]**2+self.arts.xyz.cell[0][1]**2+self.arts.xyz.cell[0][2]**2)
                b = np.sqrt(self.arts.xyz.cell[1][0]**2+self.arts.xyz.cell[1][1]**2+self.arts.xyz.cell[1][2]**2)
                c = np.sqrt(self.arts.xyz.cell[2][0]**2+self.arts.xyz.cell[2][1]**2+self.arts.xyz.cell[2][2]**2)

                fout.write("a_in=%f\n" % a)
                fout.write("b_in=%f\n" % b)
                fout.write("c_in=%f\n" % c)

                fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])

                range_a_start = range_a[0] + i_batch_a * self.batch_a * range_a[2]
                range_a_end = range_a[0] + (i_batch_a+1) * self.batch_a * range_a[2] - range_a[2] / 2
                # - range_a[2] / 2, so that the last value is ignored which is actually the begining of next batch
                if range_a_end  > range_a[1]:
                    range_a_end = range_a[1]

                fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                fout.write("do\n")
                fout.write("  cp relax.in.template relax-${a}.in\n")
                fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${a} / ${c_in}; print result\" | bc`)\n")
                fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${a} / ${c_in}; print result\" | bc`)\n")                                        
                fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${a} / ${c_in}; print result\" | bc`)\n")                                          
                fout.write("  cat >> relax-${a}.in <<EOF\n")
                fout.write("\n")
                fout.write("CELL_PARAMETERS angstrom\n")
                fout.write("${vec11} ${vec12} ${vec13}\n")
                fout.write("${vec21} ${vec22} ${vec23}\n")
                fout.write("${vec31} ${vec32} ${vec33}\n")
                fout.write("EOF\n")
                fout.write("  yhrun $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
                fout.write("done\n")


            # gen pbs script
            with open("relax-cubic-%d.pbs" % (i_batch_a), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s-%d\n" % (self.run_params["jobname"], i_batch_a))
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])            
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

                a = np.sqrt(self.arts.xyz.cell[0][0]**2+self.arts.xyz.cell[0][1]**2+self.arts.xyz.cell[0][2]**2)
                b = np.sqrt(self.arts.xyz.cell[1][0]**2+self.arts.xyz.cell[1][1]**2+self.arts.xyz.cell[1][2]**2)
                c = np.sqrt(self.arts.xyz.cell[2][0]**2+self.arts.xyz.cell[2][1]**2+self.arts.xyz.cell[2][2]**2)

                fout.write("a_in=%f\n" % a)
                fout.write("b_in=%f\n" % b)
                fout.write("c_in=%f\n" % c)

                fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])
                
                range_a_start = range_a[0] + i_batch_a * self.batch_a * range_a[2]
                range_a_end = range_a[0] + (i_batch_a+1) * self.batch_a * range_a[2] - range_a[2] / 2
                # - range_a[2] / 2, so that the last value is ignored which is actually the begining of next batch
                if range_a_end  > range_a[1]:
                    range_a_end = range_a[1]

                fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                fout.write("do\n")
                fout.write("  cp relax.in.template relax-${a}.in\n")
                fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${a} / ${c_in}; print result\" | bc`)\n")
                fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${a} / ${c_in}; print result\" | bc`)\n")                                        
                fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${a} / ${c_in}; print result\" | bc`)\n")                              
                fout.write("  cat >> relax-${a}.in <<EOF\n")
                fout.write("\n")
                fout.write("CELL_PARAMETERS angstrom\n")
                fout.write("${vec11} ${vec12} ${vec13}\n")
                fout.write("${vec21} ${vec22} ${vec23}\n")
                fout.write("${vec31} ${vec32} ${vec33}\n")
                fout.write("EOF\n")
                fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
                fout.write("done\n")

            # gen local bash script
            with open("relax-cubic-%d.sh" % (i_batch_a), 'w') as fout:
                fout.write("#!/bin/bash\n")

                a = np.sqrt(self.arts.xyz.cell[0][0]**2+self.arts.xyz.cell[0][1]**2+self.arts.xyz.cell[0][2]**2)
                b = np.sqrt(self.arts.xyz.cell[1][0]**2+self.arts.xyz.cell[1][1]**2+self.arts.xyz.cell[1][2]**2)
                c = np.sqrt(self.arts.xyz.cell[2][0]**2+self.arts.xyz.cell[2][1]**2+self.arts.xyz.cell[2][2]**2)

                fout.write("a_in=%f\n" % a)
                fout.write("b_in=%f\n" % b)
                fout.write("c_in=%f\n" % c)

                fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])

                range_a_start = range_a[0] + i_batch_a * self.batch_a * range_a[2]
                range_a_end = range_a[0] + (i_batch_a+1) * self.batch_a * range_a[2] - range_a[2] / 2
                # - range_a[2] / 2, so that the last value is ignored which is actually the begining of next batch
                if range_a_end  > range_a[1]:
                    range_a_end = range_a[1]

                fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                fout.write("do\n")
                fout.write("  cp relax.in.template relax-${a}.in\n")
                fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${a} / ${c_in}; print result\" | bc`)\n")
                fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${a} / ${c_in}; print result\" | bc`)\n")                                        
                fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${a} / ${c_in}; print result\" | bc`)\n")                              
                fout.write("  cat >> relax-${a}.in <<EOF\n")
                fout.write("\n")
                fout.write("CELL_PARAMETERS angstrom\n")
                fout.write("${vec11} ${vec12} ${vec13}\n")
                fout.write("${vec21} ${vec22} ${vec23}\n")
                fout.write("${vec31} ${vec32} ${vec33}\n")
                fout.write("EOF\n")
                fout.write("  %s $PMF_PWX < relax-${a}.in | tee relax-${a}.out\n" % self.run_params["mpi"])
                fout.write("done\n")

            # gen cdcloud script
            with open("relax-cubic-%d.slurm_cd" % (i_batch_a), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s-%d\n" % (self.run_params["jobname"], i_batch_a))
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                fout.write("#\n")
                fout.write("export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so\n")           
                #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

                a = np.sqrt(self.arts.xyz.cell[0][0]**2+self.arts.xyz.cell[0][1]**2+self.arts.xyz.cell[0][2]**2)
                b = np.sqrt(self.arts.xyz.cell[1][0]**2+self.arts.xyz.cell[1][1]**2+self.arts.xyz.cell[1][2]**2)
                c = np.sqrt(self.arts.xyz.cell[2][0]**2+self.arts.xyz.cell[2][1]**2+self.arts.xyz.cell[2][2]**2)

                fout.write("a_in=%f\n" % a)
                fout.write("b_in=%f\n" % b)
                fout.write("c_in=%f\n" % c)

                fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])

                range_a_start = range_a[0] + i_batch_a * self.batch_a * range_a[2]
                range_a_end = range_a[0] + (i_batch_a+1) * self.batch_a * range_a[2] - range_a[2] / 2
                # - range_a[2] / 2, so that the last value is ignored which is actually the begining of next batch
                if range_a_end  > range_a[1]:
                    range_a_end = range_a[1]

                fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                fout.write("do\n")
                fout.write("  cp relax.in.template relax-${a}.in\n")
                fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${a} / ${c_in}; print result\" | bc`)\n")
                fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${a} / ${c_in}; print result\" | bc`)\n")                                        
                fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${a} / ${c_in}; print result\" | bc`)\n")                                          
                fout.write("  cat >> relax-${a}.in <<EOF\n")
                fout.write("\n")
                fout.write("CELL_PARAMETERS angstrom\n")
                fout.write("${vec11} ${vec12} ${vec13}\n")
                fout.write("${vec21} ${vec22} ${vec23}\n")
                fout.write("${vec31} ${vec32} ${vec33}\n")
                fout.write("EOF\n")
                fout.write("  srun --mpi=pmix_v3 $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
                fout.write("done\n")

        # generate result analysis script
        os.system("mkdir -p post-processing")

        with open("post-processing/get_energy.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("cat > energy-latconst.data <<EOF\n")
            fout.write("# format: a energy(Ry)\n")
            fout.write("EOF\n")
            fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a[0], range_a[2], a+range_a[1]))
            fout.write("do\n")
            fout.write("  energy=`cat ../relax-${a}.out | grep '!    total energy' | tail -1`\n")
            fout.write("  cat >> energy-latconst.data <<EOF\n")
            fout.write("${a} ${energy:32:-2}\n")
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
            for i_batch_a in range(n_batch_a):
                os.system("bash relax-cubic-%d.sh" % i_batch_a)
            os.chdir("../")
        for i_batch_a in range(n_batch_a):
            server_handle(auto=auto, directory=directory, jobfilebase="relax-cubic-%d" % i_batch_a, server=self.run_params["server"])

    def hexagonal(self, directory="tmp-qe-hexagonal", runopt="gen", auto=0, range_a=[-0.1, 0.101, 0.01], range_c=[-0.1, 0.101, 0.01]):
        """
        """
        na = len(np.arange(range_a[0], range_a[1], range_a[2]))
        nc = len(np.arange(range_c[0], range_c[1], range_c[2]))

        if self.batch_a == None:
            # namely all in one batch
            self.batch_a = na
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

        if nc % self.batch_c == 0:
            n_batch_c = int(nc / self.batch_c)
        else:
            n_batch_c = int(nc / self.batch_c) + 1        
        #

        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)

        shutil.copyfile(self.arts.xyz.file, os.path.join(directory, os.path.basename(self.arts.xyz.file)))
        #all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
        all_file = os.listdir()
        for element in self.arts.xyz.specie_labels:
            for item in all_file:
                #if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                #if item.split(".")[0].lower() == element.lower() or item.split("_")[0].lower() == element.lower():
                if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE) or re.match("(%s)(_*)(upf)" % element, item, re.IGNORECASE):    
                    shutil.copyfile(item, os.path.join(directory, item))
                    break
        self.arts.pseudo.dir = os.path.abspath(directory)
        self.control.set_params({"pseudo_dir": os.path.abspath(directory)})
        #
        os.chdir(directory)

        with open("relax.in.template", 'w') as fout:
            self.control.to_in(fout)
            self.system.to_in(fout)
            self.electrons.to_in(fout)
            self.ions.to_in(fout)

            coordtype = "crystal" # use crystal here so we could only change cell when opt cell
            fout.write("ATOMIC_SPECIES\n")
            all_file = os.listdir(self.arts.pseudo.dir)
            for element in self.arts.xyz.specie_labels:
                for item in all_file:
                    if re.match("(%s)(.*)(upf)" % (element), item, re.IGNORECASE):
                        fout.write("%s %f %s\n" % (element, base.element[element].mass, item))
                        break
            fout.write("\n")
            if coordtype == "angstrom":
                fout.write("ATOMIC_POSITIONS angstrom\n")
                if self.arts.ifstatic == True:
                    for atom in self.arts.xyz.atoms:
                        fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
                elif self.arts.ifstatic == False:
                    for atom in self.arts.xyz.atoms:
                        fout.write("%s\t%.9f\t%.9f\t%.9f" % (atom.name, atom.x, atom.y, atom.z))
                        for fix in atom.fix:
                            if fix == True:
                                fout.write("\t0")
                            elif fix == False:
                                fout.write("\t1")
                        fout.write("\n")
                else:
                    print("===============================================\n")
                    print("warning: qe.base.arts.to_in():\n")
                    print("arts.ifstatic could only be True or False\n")
                    sys.exit(1)
                fout.write("\n")
            elif coordtype == "crystal":
                # crystal namely fractional coordinate can be convert from cartesian coordinates
                # the conversion process is like transformation of presentation in quantum mechanics
                # the convmat is bulid to do the conversion
                #latcell = np.array(self.xyz.cell)
                #latcell = latcell.reshape(3, 3)
                latcell = np.array(self.arts.xyz.cell)
                convmat = np.linalg.inv(latcell.T)
                crystal_coord = np.zeros([self.arts.xyz.natom, 3])
                for i in range(self.arts.xyz.natom):
                    crystal_coord[i] = convmat.dot(np.array([self.arts.xyz.atoms[i].x, self.arts.xyz.atoms[i].y, self.arts.xyz.atoms[i].z]))
                #
                fout.write("ATOMIC_POSITIONS crystal\n")
                if self.arts.ifstatic == True:
                    for k in range(self.arts.xyz.natom):
                        fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (self.arts.xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))
                elif self.arts.ifstatic == False:
                    for k in range(self.arts.xyz.natom):
                        fout.write("%s\t%.9f\t%.9f\t%.9f" % (self.arts.xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))
                        for fix in self.arts.xyz.atoms[k].fix:
                            if fix == True:
                                fout.write("\t0")
                            elif fix == False:
                                fout.write("\t1")
                        fout.write("\n")
                else:
                    print("===============================================\n")
                    print("warning: qe.base.arts.to_in():\n")
                    print("arts.ifstatic could only be True or False\n")
                    sys.exit(1)
                fout.write("\n")
            # end crystal type ATOMIC_POSITIONS

            # writing KPOINTS to the fout
            self.arts.write_kpoints(fout)
            # =========================
            #
            # writing forces act on atoms
            if self.arts.atomic_forces_status == True:
                self.arts.write_atomic_forces(fout)
            # =========================

        for i_batch_a in range(n_batch_a):
            for i_batch_c in range(n_batch_c):
                # gen llhpc script
                with open("relax-hexagonal-%d-%d.slurm" % (i_batch_a, i_batch_c), 'w') as fout:
                    fout.write("#!/bin/bash\n")
                    fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                    fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                    fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                    fout.write("#SBATCH -J %s-%d-%d\n" % (self.run_params["jobname"], i_batch_a, i_batch_c))
                    fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                    fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                    #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

                    a = np.sqrt(self.arts.xyz.cell[0][0]**2+self.arts.xyz.cell[0][1]**2+self.arts.xyz.cell[0][2]**2)
                    b = np.sqrt(self.arts.xyz.cell[1][0]**2+self.arts.xyz.cell[1][1]**2+self.arts.xyz.cell[1][2]**2)
                    c = np.sqrt(self.arts.xyz.cell[2][0]**2+self.arts.xyz.cell[2][1]**2+self.arts.xyz.cell[2][2]**2)

                    fout.write("a_in=%f\n" % a)
                    fout.write("b_in=%f\n" % b)
                    fout.write("c_in=%f\n" % c)

                    fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                    fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                    fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                    fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                    fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                    fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                    fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                    fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                    fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])

                    range_a_start = range_a[0] + i_batch_a * self.batch_a * range_a[2]
                    range_a_end = range_a[0] + (i_batch_a+1) * self.batch_a * range_a[2] - range_a[2] / 2
                    # - range_a[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_a_end  > range_a[1]:
                        range_a_end = range_a[1]

                    range_c_start = range_c[0] + i_batch_c * self.batch_c * range_c[2]
                    range_c_end = range_c[0] + (i_batch_c+1) * self.batch_c * range_c[2] - range_c[2] / 2
                    # - range_c[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_c_end  > range_c[1]:
                        range_c_end = range_c[1]

                    if na >= 2:
                        # a is optimized                        
                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        if nc >= 2:
                            # optimize both a and c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                            # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                            fout.write("  cp relax.in.template relax-${a}-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")                                      
                            fout.write("  cat >> relax-${a}-${c}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  yhrun $PMF_PWX < relax-${a}-${c}.in > relax-${a}-${c}.out\n")
                            fout.write("done\n")
                        else:
                            # only optimize a
                            fout.write("  cp relax.in.template relax-${a}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; print result\" | bc`)\n")                                      
                            fout.write("  cat >> relax-${a}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  yhrun $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
                        fout.write("done\n")
                    else:
                       # a is not optimized
                        if nc >= 2:
                            # only optimize c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            fout.write("  cp relax.in.template relax-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")                                      
                            fout.write("  cat >> relax-${c}.in<<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  yhrun $PMF_PWX < relax-${c}.in > relax-${c}.out\n")
                            fout.write("done\n")
                        else:
                            # neither a or c is optimized
                            pass

                # gen pbs script
                with open("relax-hexagonal-%d-%d.pbs" % (i_batch_a, i_batch_c), 'w') as fout:
                    fout.write("#!/bin/bash\n")
                    fout.write("#PBS -N %s-%d-%d\n" % (self.run_params["jobname"], i_batch_a, i_batch_c))
                    fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                    if "queue" in self.run_params and self.run_params["queue"] != None:
                        fout.write("#PBS -q %s\n" %self.run_params["queue"])            
                    fout.write("\n")
                    fout.write("cd $PBS_O_WORKDIR\n")
                    fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                    #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

                    a = np.sqrt(self.arts.xyz.cell[0][0]**2+self.arts.xyz.cell[0][1]**2+self.arts.xyz.cell[0][2]**2)
                    b = np.sqrt(self.arts.xyz.cell[1][0]**2+self.arts.xyz.cell[1][1]**2+self.arts.xyz.cell[1][2]**2)
                    c = np.sqrt(self.arts.xyz.cell[2][0]**2+self.arts.xyz.cell[2][1]**2+self.arts.xyz.cell[2][2]**2)

                    fout.write("a_in=%f\n" % a)
                    fout.write("b_in=%f\n" % b)
                    fout.write("c_in=%f\n" % c)

                    fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                    fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                    fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                    fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                    fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                    fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                    fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                    fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                    fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])

                    range_a_start = range_a[0] + i_batch_a * self.batch_a * range_a[2]
                    range_a_end = range_a[0] + (i_batch_a+1) * self.batch_a * range_a[2] - range_a[2] / 2
                    # - range_a[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_a_end  > range_a[1]:
                        range_a_end = range_a[1]

                    range_c_start = range_c[0] + i_batch_c * self.batch_c * range_c[2]
                    range_c_end = range_c[0] + (i_batch_c+1) * self.batch_c * range_c[2] - range_c[2] / 2
                    # - range_c[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_c_end  > range_c[1]:
                        range_c_end = range_c[1]

                    if na >= 2:
                        # a is optimized
                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        if nc >= 2:
                            # optimize both a and c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                            # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                            fout.write("  cp relax.in.template relax-${a}-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")                                      
                            fout.write("  cat >> relax-${a}-${c}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < relax-${a}-${c}.in > relax-${a}-${c}.out\n")
                            fout.write("done\n")
                        else:
                            # only optimize a
                            fout.write("  cp relax.in.template relax-${a}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; print result\" | bc`)\n")                                      
                            fout.write("  cat >> relax-${a}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
                        fout.write("done\n")
                    else:
                        # a is not optimized
                        if nc >= 2:
                            # only optimize c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], range_c_end))
                            fout.write("do\n")
                            fout.write("  cp relax.in.template relax-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")                                      
                            fout.write("  cat >> relax-${c}.in<<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < relax-${c}.in > relax-${c}.out\n")
                            fout.write("done\n")
                        else:
                            # neither a or c is optimized
                            pass

                # gen local bash script
                with open("relax-hexagonal-%d-%d.sh" % (i_batch_a, i_batch_c), 'w') as fout:
                    fout.write("#!/bin/bash\n")

                    a = np.sqrt(self.arts.xyz.cell[0][0]**2+self.arts.xyz.cell[0][1]**2+self.arts.xyz.cell[0][2]**2)
                    b = np.sqrt(self.arts.xyz.cell[1][0]**2+self.arts.xyz.cell[1][1]**2+self.arts.xyz.cell[1][2]**2)
                    c = np.sqrt(self.arts.xyz.cell[2][0]**2+self.arts.xyz.cell[2][1]**2+self.arts.xyz.cell[2][2]**2)

                    fout.write("a_in=%f\n" % a)
                    fout.write("b_in=%f\n" % b)
                    fout.write("c_in=%f\n" % c)

                    fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                    fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                    fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                    fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                    fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                    fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                    fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                    fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                    fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])

                    range_a_start = range_a[0] + i_batch_a * self.batch_a * range_a[2]
                    range_a_end = range_a[0] + (i_batch_a+1) * self.batch_a * range_a[2] - range_a[2] / 2
                    # - range_a[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_a_end  > range_a[1]:
                        range_a_end = range_a[1]

                    range_c_start = range_c[0] + i_batch_c * self.batch_c * range_c[2]
                    range_c_end = range_c[0] + (i_batch_c+1) * self.batch_c * range_c[2] - range_c[2] / 2
                    # - range_c[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_c_end  > range_c[1]:
                        range_c_end = range_c[1]

                    if na >= 2:
                        # a is optimized
                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        if nc >= 2:
                            # optimize both a and c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                            # set as '0.123', namely the ommited 0 by bs by default is not ommited now!                    
                            fout.write("  cp relax.in.template relax-${a}-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")                                  
                            fout.write("  cat >> relax-${a}-${c}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  %s $PMF_PWX < relax-${a}-${c}.in | tee relax-${a}-${c}.out\n" % self.run_params["mpi"])
                            fout.write("done\n")
                        else:
                            # only optimize a
                            fout.write("  cp relax.in.template relax-${a}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; print result\" | bc`)\n")                                  
                            fout.write("  cat >> relax-${a}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  %s $PMF_PWX < relax-${a}.in | tee relax-${a}.out\n" % self.run_params["mpi"])
                        fout.write("done\n")
                    else:
                        # a is not optimized
                        if nc >= 2:
                            # only optimize c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            fout.write("  cp relax.in.template relax-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")              
                            fout.write("  cat >> relax-${c}.in<<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  %s $PMF_PWX < relax-${c}.in | tee relax-${c}.out\n" % self.run_params["mpi"])
                            fout.write("done\n")
                        else:
                            # neither a or c is optimized
                            pass

                # gen cdcloud script
                with open("relax-hexagonal-%d-%d.slurm_cd" % (i_batch_a, i_batch_c), 'w') as fout:
                    fout.write("#!/bin/bash\n")
                    fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                    fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                    fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                    fout.write("#SBATCH -J %s-%d-%d\n" % (self.run_params["jobname"], i_batch_a, i_batch_c))
                    fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                    fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                    fout.write("#\n")
                    fout.write("export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so\n")
                    #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

                    a = np.sqrt(self.arts.xyz.cell[0][0]**2+self.arts.xyz.cell[0][1]**2+self.arts.xyz.cell[0][2]**2)
                    b = np.sqrt(self.arts.xyz.cell[1][0]**2+self.arts.xyz.cell[1][1]**2+self.arts.xyz.cell[1][2]**2)
                    c = np.sqrt(self.arts.xyz.cell[2][0]**2+self.arts.xyz.cell[2][1]**2+self.arts.xyz.cell[2][2]**2)

                    fout.write("a_in=%f\n" % a)
                    fout.write("b_in=%f\n" % b)
                    fout.write("c_in=%f\n" % c)

                    fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                    fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                    fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                    fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                    fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                    fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                    fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                    fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                    fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])

                    range_a_start = range_a[0] + i_batch_a * self.batch_a * range_a[2]
                    range_a_end = range_a[0] + (i_batch_a+1) * self.batch_a * range_a[2] - range_a[2] / 2
                    # - range_a[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_a_end  > range_a[1]:
                        range_a_end = range_a[1]

                    range_c_start = range_c[0] + i_batch_c * self.batch_c * range_c[2]
                    range_c_end = range_c[0] + (i_batch_c+1) * self.batch_c * range_c[2] - range_c[2] / 2
                    # - range_c[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_c_end  > range_c[1]:
                        range_c_end = range_c[1]

                    if na >= 2:
                        # a is optimized                        
                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        if nc >= 2:
                            # optimize both a and c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                            # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                            fout.write("  cp relax.in.template relax-${a}-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")                                      
                            fout.write("  cat >> relax-${a}-${c}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  srun --mpi=pmix_v3 $PMF_PWX < relax-${a}-${c}.in > relax-${a}-${c}.out\n")
                            fout.write("done\n")
                        else:
                            # only optimize a
                            fout.write("  cp relax.in.template relax-${a}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; print result\" | bc`)\n")                                      
                            fout.write("  cat >> relax-${a}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  srun --mpi=pmix_v3 $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
                        fout.write("done\n")
                    else:
                       # a is not optimized
                        if nc >= 2:
                            # only optimize c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            fout.write("  cp relax.in.template relax-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")                                      
                            fout.write("  cat >> relax-${c}.in<<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  srun --mpi=pmix_v3 $PMF_PWX < relax-${c}.in > relax-${c}.out\n")
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
                fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a[0], range_a[2], a+range_a[1]))
                fout.write("do\n")
                if nc >= 2:
                    # both a and c are optimized
                    fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c[0], range_c[2], c+range_c[1]))
                    fout.write("do\n")
                    fout.write("  energy=`cat ../relax-${a}-${c}.out | grep '!    total energy' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${c} ${energy:32:-2}\n")
                    fout.write("EOF\n")
                    fout.write("done\n")
                    fout.write("doen\n")
                    fout.write("cat > energy-latconst.gp<<EOF\n")
                    fout.write("set term gif\n")
                    fout.write("set output 'energy-latconst.gif'\n")
                    fout.write("set title 'Energy Latconst'\n")
                    fout.write("set xlabel 'latconst(a)'\n")
                    fout.write("set ylabel 'latconst(c)'\n")
                    fout.write("set zlabel 'Energy'\n")
                    fout.write("splot 'energy-latconst.data'\n")
                    fout.write("EOF\n")
                else:
                    fout.write("  energy=`cat ../relax-${a}.out | grep '!    total energy' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${energy:32:-2}\n")
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
            else:
                # a is not optimized
                if nc >= 2:
                    # only c is optimized
                    fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c[0], range_c[2], c+range_c[1]))
                    fout.write("do\n")
                    fout.write("  energy=`cat ../relax-${c}.out | grep '!    total energy' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${c} ${energy:32:-2}\n")
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
                else:
                    # neither a nor c is optimized
                    pass
            fout.write("gnuplot energy-latconst.gp\n")                
        #os.system("cd post-processing; bash get_energy.sh; cd ../")
        os.chdir("../")
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for i_batch_a in range(n_batch_a):
                for i_batch_c in range(n_batch_c):
                    os.system("bash relax-hexagonal-%d-%d.sh" % (i_batch_a, i_batch_c))
            os.chdir("../")
        for i_batch_a in range(n_batch_a):
            for i_batch_c in range(n_batch_c):
                server_handle(auto=auto, directory=directory, jobfilebase="relax-hexagonal-%d-%d" % (i_batch_a, i_batch_c), server=self.run_params["server"])

    def tetragonal(self, directory="tmp-qe-relax-tetragonal", runopt="gen",  auto=0, range_a=[-0.1, 0.101, 0.01], range_c=[-0.1, 0.101, 0.01]):
        """
        """
        na = len(np.arange(range_a[0], range_a[1], range_a[2]))
        nc = len(np.arange(range_c[0], range_c[1], range_c[2]))

        if self.batch_a == None:
            # namely all in one batch
            self.batch_a = na
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

        if nc % self.batch_c == 0:
            n_batch_c = int(nc / self.batch_c)
        else:
            n_batch_c = int(nc / self.batch_c) + 1        
        #

        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)

        shutil.copyfile(self.arts.xyz.file, os.path.join(directory, os.path.basename(self.arts.xyz.file)))
        #all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
        all_file = os.listdir()
        for element in self.arts.xyz.specie_labels:
            for item in all_file:
                #if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                #if item.split(".")[0].lower() == element.lower() or item.split("_")[0].lower() == element.lower():
                if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE) or re.match("(%s)(_*)(upf)" % element, item, re.IGNORECASE):    
                    shutil.copyfile(item, os.path.join(directory, item))
                    break
        self.arts.pseudo.dir = os.path.abspath(directory)
        self.control.set_params({"pseudo_dir": os.path.abspath(directory)})
        #
        os.chdir(directory)

        with open("relax.in.template", 'w') as fout:
            self.control.to_in(fout)
            self.system.to_in(fout)
            self.electrons.to_in(fout)
            self.ions.to_in(fout)

            coordtype = "crystal" # use crystal here so we could only change cell when opt cell
            fout.write("ATOMIC_SPECIES\n")
            all_file = os.listdir(self.arts.pseudo.dir)
            for element in self.arts.xyz.specie_labels:
                for item in all_file:
                    if re.match("(%s)(.*)(upf)" % (element), item, re.IGNORECASE):
                        fout.write("%s %f %s\n" % (element, base.element[element].mass, item))
                        break
            fout.write("\n")
            if coordtype == "angstrom":
                fout.write("ATOMIC_POSITIONS angstrom\n")
                if self.arts.ifstatic == True:
                    for atom in self.arts.xyz.atoms:
                        fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
                elif self.arts.ifstatic == False:
                    for atom in self.arts.xyz.atoms:
                        fout.write("%s\t%.9f\t%.9f\t%.9f" % (atom.name, atom.x, atom.y, atom.z))
                        for fix in atom.fix:
                            if fix == True:
                                fout.write("\t0")
                            elif fix == False:
                                fout.write("\t1")
                        fout.write("\n")
                else:
                    print("===============================================\n")
                    print("warning: qe.base.arts.to_in():\n")
                    print("arts.ifstatic could only be True or False\n")
                    sys.exit(1)
                fout.write("\n")
            elif coordtype == "crystal":
                # crystal namely fractional coordinate can be convert from cartesian coordinates
                # the conversion process is like transformation of presentation in quantum mechanics
                # the convmat is bulid to do the conversion
                #latcell = np.array(self.xyz.cell)
                #latcell = latcell.reshape(3, 3)
                latcell = np.array(self.arts.xyz.cell)
                convmat = np.linalg.inv(latcell.T)
                crystal_coord = np.zeros([self.arts.xyz.natom, 3])
                for i in range(self.arts.xyz.natom):
                    crystal_coord[i] = convmat.dot(np.array([self.arts.xyz.atoms[i].x, self.arts.xyz.atoms[i].y, self.arts.xyz.atoms[i].z]))
                #
                fout.write("ATOMIC_POSITIONS crystal\n")
                if self.arts.ifstatic == True:
                    for k in range(self.arts.xyz.natom):
                        fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (self.arts.xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))
                elif self.arts.ifstatic == False:
                    for k in range(self.arts.xyz.natom):
                        fout.write("%s\t%.9f\t%.9f\t%.9f" % (self.arts.xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))
                        for fix in self.arts.xyz.atoms[k].fix:
                            if fix == True:
                                fout.write("\t0")
                            elif fix == False:
                                fout.write("\t1")
                        fout.write("\n")
                else:
                    print("===============================================\n")
                    print("warning: qe.base.arts.to_in():\n")
                    print("arts.ifstatic could only be True or False\n")
                    sys.exit(1)
                fout.write("\n")
            # end crystal type ATOMIC_POSITIONS

            # writing KPOINTS to the fout
            self.arts.write_kpoints(fout)
            # =========================
            #
            # writing forces act on atoms
            if self.arts.atomic_forces_status == True:
                self.arts.write_atomic_forces(fout)
            # =========================
        for i_batch_a in range(n_batch_a):
            for i_batch_c in range(n_batch_c):
                # gen llhpc script
                with open("relax-tetragonal-%d-%d.slurm" % (i_batch_a, i_batch_c), 'w') as fout:
                    fout.write("#!/bin/bash\n")
                    fout.write("#!/bin/bash\n")
                    fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                    fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                    fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                    fout.write("#SBATCH -J %s-%d-%d\n" % (self.run_params["jobname"], i_batch_a, i_batch_c))
                    fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                    fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                    #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

                    a = np.sqrt(self.arts.xyz.cell[0][0]**2+self.arts.xyz.cell[0][1]**2+self.arts.xyz.cell[0][2]**2)
                    b = np.sqrt(self.arts.xyz.cell[1][0]**2+self.arts.xyz.cell[1][1]**2+self.arts.xyz.cell[1][2]**2)
                    c = np.sqrt(self.arts.xyz.cell[2][0]**2+self.arts.xyz.cell[2][1]**2+self.arts.xyz.cell[2][2]**2)

                    fout.write("a_in=%f\n" % a)
                    fout.write("b_in=%f\n" % b)
                    fout.write("c_in=%f\n" % c)

                    fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                    fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                    fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                    fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                    fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                    fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                    fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                    fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                    fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])

                    range_a_start = range_a[0] + i_batch_a * self.batch_a * range_a[2]
                    range_a_end = range_a[0] + (i_batch_a+1) * self.batch_a * range_a[2] - range_a[2] / 2
                    # - range_a[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_a_end  > range_a[1]:
                        range_a_end = range_a[1]

                    range_c_start = range_c[0] + i_batch_c * self.batch_c * range_c[2]
                    range_c_end = range_c[0] + (i_batch_c+1) * self.batch_c * range_c[2] - range_c[2] / 2
                    # - range_c[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_c_end  > range_c[1]:
                        range_c_end = range_c[1]

                    if na >= 2:
                        # a is optimized
                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        if nc >= 2:
                            # optimize both a and c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            fout.write("    cp relax.in.template relax-${a}-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")                                  
                            fout.write("    cat >> relax-${a}-${c}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("    yhrun $PMF_PWX < relax-${a}-${c}.in > relax-${a}-${c}.out\n")
                            fout.write("  done\n")
                        else:
                            # only optimize a
                            fout.write("    cp relax.in.template relax-${a}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; print result\" | bc`)\n")                                  
                            fout.write("    cat >> relax-${a}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("    yhrun $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
                        fout.write("done\n")
                    else:
                        # a is not optimized
                        if nc >= 2:
                            # only optimize c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            fout.write("  cp relax.in.template relax-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")              
                            fout.write("  cat >> relax-${c}.in<<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  yhrun $PMF_PWX < relax-${c}.in > relax-${c}.out\n")
                            fout.write("done\n")
                        else:
                            # neither a or c is optimized
                            pass


                # gen pbs script
                with open("relax-tetragonal-%d-%d.pbs" % (i_batch_a, i_batch_c), 'w') as fout:
                    fout.write("#!/bin/bash\n")
                    fout.write("#PBS -N %s-%d-%d\n" % (self.run_params["jobname"], i_batch_a, i_batch_c))
                    fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                    if "queue" in self.run_params and self.run_params["queue"] != None:
                        fout.write("#PBS -q %s\n" %self.run_params["queue"])            
                    fout.write("\n")
                    fout.write("cd $PBS_O_WORKDIR\n")
                    fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                    #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

                    a = np.sqrt(self.arts.xyz.cell[0][0]**2+self.arts.xyz.cell[0][1]**2+self.arts.xyz.cell[0][2]**2)
                    b = np.sqrt(self.arts.xyz.cell[1][0]**2+self.arts.xyz.cell[1][1]**2+self.arts.xyz.cell[1][2]**2)
                    c = np.sqrt(self.arts.xyz.cell[2][0]**2+self.arts.xyz.cell[2][1]**2+self.arts.xyz.cell[2][2]**2)

                    fout.write("a_in=%f\n" % a)
                    fout.write("b_in=%f\n" % b)
                    fout.write("c_in=%f\n" % c)

                    fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                    fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                    fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                    fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                    fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                    fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                    fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                    fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                    fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])

                    range_a_start = range_a[0] + i_batch_a * self.batch_a * range_a[2]
                    range_a_end = range_a[0] + (i_batch_a+1) * self.batch_a * range_a[2] - range_a[2] / 2
                    # - range_a[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_a_end  > range_a[1]:
                        range_a_end = range_a[1]

                    range_c_start = range_c[0] + i_batch_c * self.batch_c * range_c[2]
                    range_c_end = range_c[0] + (i_batch_c+1) * self.batch_c * range_c[2] - range_c[2] / 2
                    # - range_c[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_c_end  > range_c[1]:
                        range_c_end = range_c[1]

                    if na >= 2:
                       # a is optimized
                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        if nc >= 2:
                            # optimize both a and c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            fout.write("    cp relax.in.template relax-${a}-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")              
                            fout.write("    cat >> relax-${a}-${c}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("    mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < relax-${a}-${c}.in > relax-${a}-${c}.out\n")
                            fout.write("  done\n")
                        else:
                            # only optimize a
                            fout.write("    cp relax.in.template relax-${a}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; print result\" | bc`)\n")              
                            fout.write("    cat >> relax-${a}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("    mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
                        fout.write("done\n")
                    else:
                        # a is not optimized
                        if nc >= 2:
                            # only optimize c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            fout.write("  cp relax.in.template relax-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")              
                            fout.write("  cat >> relax-${c}.in<<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < relax-${c}.in > relax-${c}.out\n")
                            fout.write("done\n")
                        else:
                            # neither a or c is optimized
                            pass           
                            
                # gen local bash script
                with open("relax-tetragonal-%d-%d.bash" % (i_batch_a, i_batch_c), 'w') as fout:
                    fout.write("#!/bin/bash\n")

                    a = np.sqrt(self.arts.xyz.cell[0][0]**2+self.arts.xyz.cell[0][1]**2+self.arts.xyz.cell[0][2]**2)
                    b = np.sqrt(self.arts.xyz.cell[1][0]**2+self.arts.xyz.cell[1][1]**2+self.arts.xyz.cell[1][2]**2)
                    c = np.sqrt(self.arts.xyz.cell[2][0]**2+self.arts.xyz.cell[2][1]**2+self.arts.xyz.cell[2][2]**2)

                    fout.write("a_in=%f\n" % a)
                    fout.write("b_in=%f\n" % b)
                    fout.write("c_in=%f\n" % c)

                    fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                    fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                    fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                    fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                    fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                    fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                    fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                    fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                    fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])

                    range_a_start = range_a[0] + i_batch_a * self.batch_a * range_a[2]
                    range_a_end = range_a[0] + (i_batch_a+1) * self.batch_a * range_a[2] - range_a[2] / 2
                    # - range_a[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_a_end  > range_a[1]:
                        range_a_end = range_a[1]

                    range_c_start = range_c[0] + i_batch_c * self.batch_c * range_c[2]
                    range_c_end = range_c[0] + (i_batch_c+1) * self.batch_c * range_c[2] - range_c[2] / 2
                    # - range_c[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_c_end  > range_c[1]:
                        range_c_end = range_c[1]

                    if na >= 2:
                        # a is optimized
                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        if nc >= 2:
                            # optimize both a and c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            fout.write("    cp relax.in.template relax-${a}-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")              
                            fout.write("    cat >> relax-${a}-${c}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("    %s $PMF_PWX < relax-${a}-${c}.in | tee relax-${a}-${c}.out\n" % self.run_params["mpi"])
                            fout.write("  done\n")
                        else:
                            # only optimize a
                            fout.write("    cp relax.in.template relax-${a}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; print result\" | bc`)\n")              
                            fout.write("    cat >> relax-${a}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("    %s $PMF_PWX < relax-${a}.in | tee relax-${a}.out\n" % self.run_params["mpi"])
                        fout.write("done\n")
                    else:
                        # a is not optimized
                        if nc >= 2:
                            # only optimize c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            fout.write("  cp relax.in.template relax-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")              
                            fout.write("  cat >> relax-${c}.in<<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  %s $PMF_PWX < relax-${c}.in | tee relax-${c}.out\n" % self.run_params["mpi"])
                            fout.write("done\n")
                        else:
                            # neither a or c is optimized
                            pass

                # gen cdcloud script
                with open("relax-tetragonal-%d-%d.slurm_cd" % (i_batch_a, i_batch_c), 'w') as fout:
                    fout.write("#!/bin/bash\n")
                    fout.write("#!/bin/bash\n")
                    fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                    fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                    fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                    fout.write("#SBATCH -J %s-%d-%d\n" % (self.run_params["jobname"], i_batch_a, i_batch_c))
                    fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                    fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                    fout.write("#\n")
                    fout.write("export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so\n")
                    #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

                    a = np.sqrt(self.arts.xyz.cell[0][0]**2+self.arts.xyz.cell[0][1]**2+self.arts.xyz.cell[0][2]**2)
                    b = np.sqrt(self.arts.xyz.cell[1][0]**2+self.arts.xyz.cell[1][1]**2+self.arts.xyz.cell[1][2]**2)
                    c = np.sqrt(self.arts.xyz.cell[2][0]**2+self.arts.xyz.cell[2][1]**2+self.arts.xyz.cell[2][2]**2)

                    fout.write("a_in=%f\n" % a)
                    fout.write("b_in=%f\n" % b)
                    fout.write("c_in=%f\n" % c)

                    fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                    fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                    fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                    fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                    fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                    fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                    fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                    fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                    fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])

                    range_a_start = range_a[0] + i_batch_a * self.batch_a * range_a[2]
                    range_a_end = range_a[0] + (i_batch_a+1) * self.batch_a * range_a[2] - range_a[2] / 2
                    # - range_a[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_a_end  > range_a[1]:
                        range_a_end = range_a[1]

                    range_c_start = range_c[0] + i_batch_c * self.batch_c * range_c[2]
                    range_c_end = range_c[0] + (i_batch_c+1) * self.batch_c * range_c[2] - range_c[2] / 2
                    # - range_c[2] / 2, so that the last value is ignored which is actually the begining of next batch
                    if range_c_end  > range_c[1]:
                        range_c_end = range_c[1]

                    if na >= 2:
                        # a is optimized
                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        if nc >= 2:
                            # optimize both a and c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            fout.write("    cp relax.in.template relax-${a}-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")                                  
                            fout.write("    cat >> relax-${a}-${c}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("    srun --mpi=pmix_v3 $PMF_PWX < relax-${a}-${c}.in > relax-${a}-${c}.out\n")
                            fout.write("  done\n")
                        else:
                            # only optimize a
                            fout.write("    cp relax.in.template relax-${a}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c_in} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c_in} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c_in} / ${c_in}; print result\" | bc`)\n")                                  
                            fout.write("    cat >> relax-${a}.in <<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("    srun --mpi=pmix_v3 $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
                        fout.write("done\n")
                    else:
                        # a is not optimized
                        if nc >= 2:
                            # only optimize c
                            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                            fout.write("do\n")
                            fout.write("  cp relax.in.template relax-${c}.in\n")
                            fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a_in} / ${a_in}; print result\" | bc`)\n")                    
                            fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${a_in} / ${a_in}; print result\" | bc`)\n")
                            fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                            fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                        
                            fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")              
                            fout.write("  cat >> relax-${c}.in<<EOF\n")
                            fout.write("\n")
                            fout.write("CELL_PARAMETERS angstrom\n")
                            fout.write("${vec11} ${vec12} ${vec13}\n")
                            fout.write("${vec21} ${vec22} ${vec23}\n")
                            fout.write("${vec31} ${vec32} ${vec33}\n")
                            fout.write("EOF\n")
                            fout.write("  srun --mpi=pmix_v3 $PMF_PWX < relax-${c}.in > relax-${c}.out\n")
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
                fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a[0], range_a[2], a+range_a[1]))
                fout.write("do\n")
                if nc >= 2:
                    # both a and c are optimized
                    fout.write("  for c in `seq -w %f %f %f`\n" % (c+range_c[0], range_c[2], c+range_c[1]))
                    fout.write("  do\n")
                    fout.write("    energy=`cat ../relax-${a}-${c}.out | grep '!    total energy' | tail -1`\n")
                    fout.write("    cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${c} ${energy:32:-2}\n")
                    fout.write("EOF\n")
                    fout.write("  done\n")
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
                else:
                    fout.write("  energy=`cat ../relax-${a}.out | grep '!    total energy' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${energy:32:-2}\n")
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
            else:
                # a is not optimized
                if nc >= 2:
                    # only c is optimized
                    fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c[0], range_c[2], c+range_c[1]))
                    fout.write("do\n")
                    fout.write("  energy=`cat ../relax-${c}.out | grep '!    total energy' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${c} ${energy:32:-2}\n")
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
                else:
                    # neither a nor c is optimized
                    pass
            fout.write("gnuplot energy-latconst.gp\n")                
        #os.system("cd post-processing; bash get_energy.sh; cd ../")
        os.chdir("../")
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for i_batch_a in range(n_batch_a):
                for i_batch_c in range(n_batch_c):            
                    os.system("bash relax-tetragonal-%d-%d.sh" % (i_batch_a, i_batch_c))
            os.chdir("../")
        for i_batch_a in range(n_batch_a):
            for i_batch_c in range(n_batch_c):            
                server_handle(auto=auto, directory=directory, jobfilebase="relax-tetragonal-%d-%d" % (i_batch_a, i_batch_c), server=self.run_params["server"])


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

        shutil.copyfile(self.arts.xyz.file, os.path.join(directory, os.path.basename(self.arts.xyz.file)))
        #all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
        all_file = os.listdir()
        for element in self.arts.xyz.specie_labels:
            for item in all_file:
                #if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                #if item.split(".")[0].lower() == element.lower() or item.split("_")[0].lower() == element.lower():
                if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE) or re.match("(%s)(_*)(upf)" % element, item, re.IGNORECASE):    
                    shutil.copyfile(item, os.path.join(directory, item))
                    break
        self.arts.pseudo.dir = os.path.abspath(directory)
        self.control.set_params({"pseudo_dir": os.path.abspath(directory)})
        #
        os.chdir(directory)

        with open("relax.in.template", 'w') as fout:
            self.control.to_in(fout)
            self.system.to_in(fout)
            self.electrons.to_in(fout)
            self.ions.to_in(fout)

            coordtype = "crystal" # use crystal here so we could only change cell when opt cell
            fout.write("ATOMIC_SPECIES\n")
            all_file = os.listdir(self.arts.pseudo.dir)
            for element in self.arts.xyz.specie_labels:
                for item in all_file:
                    #if re.match("(%s)(.*)(upf)" % (element), item, re.IGNORECASE):
                    if item.split(".")[0].lower() == element.lower() or item.split("_")[0].lower() == element.lower():
                        fout.write("%s %f %s\n" % (element, base.element[element].mass, item))
                        break
            fout.write("\n")
            if coordtype == "angstrom":
                fout.write("ATOMIC_POSITIONS angstrom\n")
                if self.arts.ifstatic == True:
                    for atom in self.arts.xyz.atoms:
                        fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
                elif self.arts.ifstatic == False:
                    for atom in self.arts.xyz.atoms:
                        fout.write("%s\t%.9f\t%.9f\t%.9f" % (atom.name, atom.x, atom.y, atom.z))
                        for fix in atom.fix:
                            if fix == True:
                                fout.write("\t0")
                            elif fix == False:
                                fout.write("\t1")
                        fout.write("\n")
                else:
                    print("===============================================\n")
                    print("warning: qe.base.arts.to_in():\n")
                    print("arts.ifstatic could only be True or False\n")
                    sys.exit(1)
                fout.write("\n")
            elif coordtype == "crystal":
                # crystal namely fractional coordinate can be convert from cartesian coordinates
                # the conversion process is like transformation of presentation in quantum mechanics
                # the convmat is bulid to do the conversion
                #latcell = np.array(self.xyz.cell)
                #latcell = latcell.reshape(3, 3)
                latcell = np.array(self.arts.xyz.cell)
                convmat = np.linalg.inv(latcell.T)
                crystal_coord = np.zeros([self.arts.xyz.natom, 3])
                for i in range(self.arts.xyz.natom):
                    crystal_coord[i] = convmat.dot(np.array([self.arts.xyz.atoms[i].x, self.arts.xyz.atoms[i].y, self.arts.xyz.atoms[i].z]))
                #
                fout.write("ATOMIC_POSITIONS crystal\n")
                if self.arts.ifstatic == True:
                    for k in range(self.arts.xyz.natom):
                        fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (self.arts.xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))
                elif self.arts.ifstatic == False:
                    for k in range(self.arts.xyz.natom):
                        fout.write("%s\t%.9f\t%.9f\t%.9f" % (self.arts.xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))
                        for fix in self.arts.xyz.atoms[k].fix:
                            if fix == True:
                                fout.write("\t0")
                            elif fix == False:
                                fout.write("\t1")
                        fout.write("\n")
                else:
                    print("===============================================\n")
                    print("warning: qe.base.arts.to_in():\n")
                    print("arts.ifstatic could only be True or False\n")
                    sys.exit(1)
                fout.write("\n")
            # end crystal type ATOMIC_POSITIONS

            # writing KPOINTS to the fout
            self.arts.write_kpoints(fout)
            # =========================
            #
            # writing forces act on atoms
            if self.arts.atomic_forces_status == True:
                self.arts.write_atomic_forces(fout)
            # =========================
            
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

                    a = np.sqrt(self.arts.xyz.cell[0][0]**2+self.arts.xyz.cell[0][1]**2+self.arts.xyz.cell[0][2]**2)
                    b = np.sqrt(self.arts.xyz.cell[1][0]**2+self.arts.xyz.cell[1][1]**2+self.arts.xyz.cell[1][2]**2)
                    c = np.sqrt(self.arts.xyz.cell[2][0]**2+self.arts.xyz.cell[2][1]**2+self.arts.xyz.cell[2][2]**2)

                    # gen llhpc script
                    with open("opt-abc-%d-%d-%d.slurm" % (i_batch_a, i_batch_b, i_batch_c), 'w') as fout:
                        fout.write("#!/bin/bash\n")
                        fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                        fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                        fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                        fout.write("#SBATCH -J %s-%d-%d-%d\n" % (self.run_params["jobname"], i_batch_a, i_batch_b, i_batch_c))
                        fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                        fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
            
                        fout.write("a_in=%f\n" % a)
                        fout.write("b_in=%f\n" % b)
                        fout.write("c_in=%f\n" % c)

                        fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                        fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                        fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                        fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                        fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                        fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                        fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                        fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                        fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])


                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        fout.write("for b in `seq -w %f %f %f`\n" % (b+range_b_start, range_b[2], b+range_b_end))
                        fout.write("do\n")
                        fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                        fout.write("do\n")
                        #fout.write("  mkdir relax-${a}-${b}-${c}\n")
                        fout.write("  cp relax.in.template relax-${a}-${b}-${c}.in\n")
                        fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                        fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                   
                        fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                        fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                        fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                       
                        fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")                         
                        fout.write("    cat >> relax-${a}-${b}-${c}.in <<EOF\n")
                        fout.write("\n")
                        fout.write("CELL_PARAMETERS angstrom\n")
                        fout.write("${vec11} ${vec12} ${vec13}\n")
                        fout.write("${vec21} ${vec22} ${vec23}\n")
                        fout.write("${vec31} ${vec32} ${vec33}\n")
                        fout.write("EOF\n")
                        fout.write("    yhrun $PMF_PWX < relax-${a}-${b}-${c}.in > relax-${a}-${b}-${c}.out\n")
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

                        fout.write("a_in=%f\n" % a)
                        fout.write("b_in=%f\n" % b)
                        fout.write("c_in=%f\n" % c)

                        fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                        fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                        fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                        fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                        fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                        fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                        fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                        fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                        fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])


                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        fout.write("for b in `seq -w %f %f %f`\n" % (b+range_b_start, range_b[2], b+range_b_end))
                        fout.write("do\n")
                        fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                        fout.write("do\n")
                        #fout.write("  mkdir relax-${a}-${b}-${c}\n")
                        fout.write("  cp relax.in.template relax-${a}-${b}-${c}.in\n")
                        fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                        fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                   
                        fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                        fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                        fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                       
                        fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")                         
                        fout.write("    cat >> relax-${a}-${b}-${c}.in <<EOF\n")
                        fout.write("\n")
                        fout.write("CELL_PARAMETERS angstrom\n")
                        fout.write("${vec11} ${vec12} ${vec13}\n")
                        fout.write("${vec21} ${vec22} ${vec23}\n")
                        fout.write("${vec31} ${vec32} ${vec33}\n")
                        fout.write("EOF\n")
                        fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < relax-${a}-${b}-${c}.in > relax-${a}-${b}-${c}.out\n")         
                        fout.write("done\n")
                        fout.write("done\n")
                        fout.write("done\n")


                    # gen local bash script
                    with open("opt-abc-%d-%d-%d.sh" % (i_batch_a, i_batch_b, i_batch_c), 'w') as fout:
                        fout.write("#!/bin/bash\n")

                        fout.write("a_in=%f\n" % a)
                        fout.write("b_in=%f\n" % b)
                        fout.write("c_in=%f\n" % c)

                        fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                        fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                        fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                        fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                        fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                        fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                        fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                        fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                        fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])

                        

                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        fout.write("for b in `seq -w %f %f %f`\n" % (b+range_b_start, range_b[2], b+range_b_end))
                        fout.write("do\n")
                        fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                        fout.write("do\n")
                        #fout.write("  mkdir relax-${a}-${b}-${c}\n")
                        fout.write("  cp relax.in.template relax-${a}-${b}-${c}.in\n")
                        fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                        fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                   
                        fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                        fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                        fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                       
                        fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")                         
                        fout.write("    cat >> relax-${a}-${b}-${c}.in <<EOF\n")
                        fout.write("\n")
                        fout.write("CELL_PARAMETERS angstrom\n")
                        fout.write("${vec11} ${vec12} ${vec13}\n")
                        fout.write("${vec21} ${vec22} ${vec23}\n")
                        fout.write("${vec31} ${vec32} ${vec33}\n")
                        fout.write("EOF\n")
                        fout.write("  %s $PMF_PWX < relax-${a}-${b}-${c}.in | tee relax-${a}-${b}-${c}.out\n" % self.run_params["mpi"])
                        fout.write("done\n")
                        fout.write("done\n")
                        fout.write("done\n")


                    # gen lsf_sz script
                    with open("opt-abc-%d-%d-%d.lsf_sz" % (i_batch_a,i_batch_b, i_batch_c), 'w') as fout:
                        fout.write("#!/bin/bash\n")
                        fout.write("APP_NAME=%s\n" % self.run_params["queue"])
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

                        fout.write("a_in=%f\n" % a)
                        fout.write("b_in=%f\n" % b)
                        fout.write("c_in=%f\n" % c)

                        fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                        fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                        fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                        fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                        fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                        fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                        fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                        fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                        fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])



                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        fout.write("for b in `seq -w %f %f %f`\n" % (b+range_b_start, range_b[2], b+range_b_end))
                        fout.write("do\n")
                        fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                        fout.write("do\n")
                        #fout.write("  mkdir relax-${a}-${b}-${c}\n")
                        fout.write("  cp relax.in.template relax-${a}-${b}-${c}.in\n")
                        fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                        fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                   
                        fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                        fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                        fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                       
                        fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")                         
                        fout.write("    cat >> relax-${a}-${b}-${c}.in <<EOF\n")
                        fout.write("\n")
                        fout.write("CELL_PARAMETERS angstrom\n")
                        fout.write("${vec11} ${vec12} ${vec13}\n")
                        fout.write("${vec21} ${vec22} ${vec23}\n")
                        fout.write("${vec31} ${vec32} ${vec33}\n")
                        fout.write("EOF\n")
                        fout.write("  %s $PMF_PWX < relax-${a}-${b}-${c}.in | tee relax-${a}-${b}-${c}.out\n" % self.run_params["mpi"])                        
                        fout.write("done\n")
                        fout.write("done\n")
                        fout.write("done\n")

                    # gen cdcloud script
                    with open("opt-abc-%d-%d-%d.slurm_cd" % (i_batch_a, i_batch_b, i_batch_c), 'w') as fout:
                        fout.write("#!/bin/bash\n")
                        fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                        fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                        fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                        fout.write("#SBATCH -J %s-%d-%d-%d\n" % (self.run_params["jobname"], i_batch_a, i_batch_b, i_batch_c))
                        fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                        fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                        fout.write("#\n")
                        fout.write("export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so\n")

                        fout.write("a_in=%f\n" % a)
                        fout.write("b_in=%f\n" % b)
                        fout.write("c_in=%f\n" % c)

                        fout.write("a1=%f\n" % self.arts.xyz.cell[0][0])
                        fout.write("a2=%f\n" % self.arts.xyz.cell[0][1])
                        fout.write("a3=%f\n" % self.arts.xyz.cell[0][2])
                        fout.write("b1=%f\n" % self.arts.xyz.cell[1][0])
                        fout.write("b2=%f\n" % self.arts.xyz.cell[1][1])
                        fout.write("b3=%f\n" % self.arts.xyz.cell[1][2])
                        fout.write("c1=%f\n" % self.arts.xyz.cell[2][0])
                        fout.write("c2=%f\n" % self.arts.xyz.cell[2][1])
                        fout.write("c3=%f\n" % self.arts.xyz.cell[2][2])


                        fout.write("for a in `seq -w %f %f %f`\n" % (a+range_a_start, range_a[2], a+range_a_end))
                        fout.write("do\n")
                        fout.write("for b in `seq -w %f %f %f`\n" % (b+range_b_start, range_b[2], b+range_b_end))
                        fout.write("do\n")
                        fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c_start, range_c[2], c+range_c_end))
                        fout.write("do\n")
                        #fout.write("  mkdir relax-${a}-${b}-${c}\n")
                        fout.write("  cp relax.in.template relax-${a}-${b}-${c}.in\n")
                        fout.write("  vec11=$(printf \"%-.6f\" `echo \"scale=6; result=${a1} * ${a} / ${a_in}; print result\" | bc`)\n")
                        fout.write("  vec12=$(printf \"%-.6f\" `echo \"scale=6; result=${a2} * ${a} / ${a_in}; print result\" | bc`)\n")                   
                        fout.write("  vec13=$(printf \"%-.6f\" `echo \"scale=6; result=${a3} * ${a} / ${a_in}; print result\" | bc`)\n")
                        fout.write("  vec21=$(printf \"%-.6f\" `echo \"scale=6; result=${b1} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec22=$(printf \"%-.6f\" `echo \"scale=6; result=${b2} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec23=$(printf \"%-.6f\" `echo \"scale=6; result=${b3} * ${b} / ${b_in}; print result\" | bc`)\n")
                        fout.write("  vec31=$(printf \"%-.6f\" `echo \"scale=6; result=${c1} * ${c} / ${c_in}; print result\" | bc`)\n")
                        fout.write("  vec32=$(printf \"%-.6f\" `echo \"scale=6; result=${c2} * ${c} / ${c_in}; print result\" | bc`)\n")                                       
                        fout.write("  vec33=$(printf \"%-.6f\" `echo \"scale=6; result=${c3} * ${c} / ${c_in}; print result\" | bc`)\n")                         
                        fout.write("    cat >> relax-${a}-${b}-${c}.in <<EOF\n")
                        fout.write("\n")
                        fout.write("CELL_PARAMETERS angstrom\n")
                        fout.write("${vec11} ${vec12} ${vec13}\n")
                        fout.write("${vec21} ${vec22} ${vec23}\n")
                        fout.write("${vec31} ${vec32} ${vec33}\n")
                        fout.write("EOF\n")
                        fout.write("    srun --mpi=pmix_v3 $PMF_PWX < relax-${a}-${b}-${c}.in > relax-${a}-${b}-${c}.out\n")
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
            fout.write("for b in `seq -w %f %f %f`\n" % (b+range_b[0], range_b[2], b+range_b[1]))
            fout.write("do\n")
            fout.write("for c in `seq -w %f %f %f`\n" % (c+range_c[0], range_c[2], c+range_c[1]))
            fout.write("do\n")
            fout.write("    energy=`cat ../relax-${a}-${b}-${c}.out | grep '!    total energy' | tail -1`\n")
            fout.write("  cat >> energy-latconst.data <<EOF\n")
            fout.write("${a} ${b} ${c} ${energy:32:-2}\n")
            fout.write("EOF\n")
            fout.write("done\n")
            fout.write("done\n")
            fout.write("done\n")
            
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
