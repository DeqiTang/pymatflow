"""
Geometric Optimization calc
"""
import os
import re
import shutil

from pymatflow.remote.server import server_handle
from pymatflow.qe.pwscf import pwscf

class opt_run(pwscf):
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
                    if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                        shutil.copyfile(item, os.path.join(directory, item))
                        break
            self.arts.pseudo.dir = os.path.abspath(directory)
            self.control.pseudo_dir = os.path.abspath(directory)
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
                    if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                        shutil.copyfile(item, os.path.join(directory, item))
                        break
            #
            self.arts.pseudo.dir = os.path.abspath(directory)
            self.control.pseudo_dir = os.path.abspath(directory)


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

    def cubic(self, directory="tmp-qe-relax-cubic", runopt="gen", auto=0):
        """
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)

        shutil.copyfile(self.arts.xyz.file, os.path.join(directory, os.path.basename(self.arts.xyz.file)))
        #all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
        all_file = os.listdir()
        for element in self.arts.xyz.specie_labels:
            for item in all_file:
                if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                    shutil.copyfile(item, os.path.join(directory, item))
                    break
        self.arts.pseudo.dir = os.path.abspath(directory)
        self.control.pseudo_dir = os.path.abspath(directory)
        #
        os.chdir(directory)

        with open("relax.in.template", 'w') as fout:
            self.control.to_in(fout)
            self.system.to_in(fout)
            self.electrons.to_in(fout)
            self.ions.to_in(fout)

            coordtype = "angstrom"
            fout.write("ATOMIC_SPECIES\n")
            upf_all = [s for s in os.listdir(self.arts.pseudo.dir) if s.split(".")[-1] == "UPF"]
            for element in self.arts.xyz.specie_labels:
                for upf in upf_all:
                    if upf.split(".")[0] == element:
                        pseudo_file =upf
                        break
                fout.write("%s %f %s\n" % (element, base.element[element].mass, pseudo_file))
                pseudo_file = None
                # after pseudo_file used, set it to None to avoid it will be used in the next element
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

        # gen llhpc script
        with open("relax-cubic.slurm", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
            fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
            fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
            fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
            fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
            fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
            #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

            a = self.arts.xyz.cell[0][0]

            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  cp relax.in.template relax-${a}.in\n")
            fout.write("  cat >> relax-${a}.in <<EOF\n")
            fout.write("\n")
            fout.write("CELL_PARAMETERS angstrom\n")
            fout.write("${a} 0.000000 0.000000\n")
            fout.write("0.000000 ${a} 0.000000\n")
            fout.write("0.000000 0.000000 ${a}\n")
            fout.write("EOF\n")
            fout.write("  yhrun $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
            fout.write("done\n")


        # gen pbs script
        with open("relax-cubic.pbs", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % self.run_params["jobname"])
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
            if "queue" in self.run_params and self.run_params["queue"] != None:
                fout.write("#PBS -q %s\n" %self.run_params["queue"])            
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

            a = self.arts.xyz.cell[0][0]

            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  cp relax.in.template relax-${a}.in\n")
            fout.write("  cat >> relax-${a}.in <<EOF\n")
            fout.write("\n")
            fout.write("CELL_PARAMETERS angstrom\n")
            fout.write("${a} 0.000000 0.000000\n")
            fout.write("0.000000 ${a} 0.000000\n")
            fout.write("0.000000 0.000000 ${a}\n")
            fout.write("EOF\n")
            fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
            fout.write("done\n")

        # gen local bash script
        with open("relax-cubic.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")

            a = self.arts.xyz.cell[0][0]

            fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
            fout.write("do\n")
            fout.write("  cp relax.in.template relax-${a}.in\n")
            fout.write("  cat >> relax-${a}.in <<EOF\n")
            fout.write("\n")
            fout.write("CELL_PARAMETERS angstrom\n")
            fout.write("${a} 0.000000 0.000000\n")
            fout.write("0.000000 ${a} 0.000000\n")
            fout.write("0.000000 0.000000 ${a}\n")
            fout.write("EOF\n")
            fout.write("  %s $PMF_PWX < relax-${a}.in | tee relax-${a}.out\n" % self.run_params["mpi"])
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
            fout.write("  energy=`cat ../relax-${a}.out | grep '!    total energy' | tail -1`\n")
            fout.write("  cat >> energy-latconst.data <<EOF\n")
            fout.write("${a} ${energy:32:-2}\n")
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
            os.system("bash relax-cubic.sh")
            os.chdir("../")
        server_handle(auto=args.auto, directory=directory, jobfilebase="relax-cubic", server=args.server)

    def hexagonal(self, directory="tmp-qe-hexagonal", runopt="gen", auto=0):
        """
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)

        shutil.copyfile(self.arts.xyz.file, os.path.join(directory, os.path.basename(self.arts.xyz.file)))
        #all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
        all_file = os.listdir()
        for element in self.arts.xyz.specie_labels:
            for item in all_file:
                if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                    shutil.copyfile(item, os.path.join(directory, item))
                    break
        self.arts.pseudo.dir = os.path.abspath(directory)
        self.control.pseudo_dir = os.path.abspath(directory)
        #
        os.chdir(directory)

        with open("relax.in.template", 'w') as fout:
            self.control.to_in(fout)
            self.system.to_in(fout)
            self.electrons.to_in(fout)
            self.ions.to_in(fout)

            coordtype = "angstrom"
            fout.write("ATOMIC_SPECIES\n")
            upf_all = [s for s in os.listdir(self.arts.pseudo.dir) if s.split(".")[-1] == "UPF"]
            for element in self.arts.xyz.specie_labels:
                for upf in upf_all:
                    if upf.split(".")[0] == element:
                        pseudo_file =upf
                        break
                fout.write("%s %f %s\n" % (element, base.element[element].mass, pseudo_file))
                pseudo_file = None
                # after pseudo_file used, set it to None to avoid it will be used in the next element
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

        # gen llhpc script
        with open("relax-hexagonal.slurm", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
            fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
            fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
            fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
            fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
            fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
            #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

            a = self.arts.xyz.cell[0][0]
            c = self.arts.xyz.cell[2][2]
            fout.write("v11=%f\n" % self.arts.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.arts.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.arts.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.arts.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.arts.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.arts.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.arts.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.arts.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.arts.xyz.cell[2][2])
            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  cp relax.in.template relax-${a}-${c}.in\n")
                    fout.write("  cat >> relax-${a}-${c}.in <<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  yhrun $PMF_PWX < relax-${a}-${c}.in > relax-${a}-${c}.out\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  cp relax.in.template relax-${a}.in\n")
                    fout.write("  cat >> relax-${a}.in <<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  yhrun $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  cp relax.in.template relax-${c}.in\n")
                    fout.write("  cat >> relax-${c}.in<<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("${v21} ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  yhrun $PMF_PWX < relax-${c}.in > relax-${c}.out\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass

        # gen pbs script
        with open("relax-hexagonal.pbs", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % self.run_params["jobname"])
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
            if "queue" in self.run_params and self.run_params["queue"] != None:
                fout.write("#PBS -q %s\n" %self.run_params["queue"])            
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

            a = self.arts.xyz.cell[0][0]
            c = self.arts.xyz.cell[2][2]
            fout.write("v11=%f\n" % self.arts.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.arts.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.arts.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.arts.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.arts.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.arts.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.arts.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.arts.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.arts.xyz.cell[2][2])
            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  cp relax.in.template relax-${a}-${c}.in\n")
                    fout.write("  cat >> relax-${a}-${c}.in <<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  yhrun $PMF_PWX < relax-${a}-${c}.in > relax-${a}-${c}.out\n")
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  cp relax.in.template relax-${a}.in\n")
                    fout.write("  cat >> relax-${a}.in <<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  yhrun $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  cp relax.in.template relax-${c}.in\n")
                    fout.write("  cat >> relax-${c}.in<<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("${v21} ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  yhrun $PMF_PWX < relax-${c}.in > relax-${c}.out\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass

        # gen local bash script
        with open("relax-hexagonal.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")

            a = self.arts.xyz.cell[0][0]
            c = self.arts.xyz.cell[2][2]
            fout.write("v11=%f\n" % self.arts.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.arts.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.arts.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.arts.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.arts.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.arts.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.arts.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.arts.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.arts.xyz.cell[2][2])
            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                    # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                    fout.write("  cp relax.in.template relax-${a}-${c}.in\n")
                    fout.write("  cat >> relax-${a}-${c}.in <<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  %s $PMF_PWX < relax-${a}-${c}.in | tee relax-${a}-${c}.out\n" % self.run_params["mpi"])
                    fout.write("done\n")
                else:
                    # only optimize a
                    fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                    fout.write("  cp relax.in.template relax-${a}.in\n")
                    fout.write("  cat >> relax-${a}.in <<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("${vec21} ${vec22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("  %s $PMF_PWX < relax-${a}.in | tee relax-${a}.out\n" % self.run_params["mpi"])
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  cp relax.in.template relax-${c}.in\n")
                    fout.write("  cat >> relax-${c}.in<<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("${v21} ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  %s $PMF_PWX < relax-${c}.in | tee relax-${c}.out\n" % self.run_params["mpi"])
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
                    fout.write("  energy=`cat ../relax-${a}-${c}.out | grep '!    total energy' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${c} ${energy:32:-2}\n")
                    fout.write("EOF\n")
                    fout.write("done\n")
                    fout.write("doen\n")
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
                    fout.write("  energy=`cat ../relax-${a}.out | grep '!    total energy' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${a} ${energy:32:-2}\n")
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
                # a is not optimized
                if nc >= 2:
                    # only c is optimized
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  energy=`cat ../relax-${c}.out | grep '!    total energy' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${c} ${energy:32:-2}\n")
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
                    # neither a nor c is optimized
                    pass
        #os.system("cd post-processing; bash get_energy.sh; cd ../")
        os.chdir("../")
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash relax-hexagonal.sh")
            os.chdir("../")
        server_handle(auto=args.auto, directory=directory, jobfilebase="relax-hexagonal", server=args.server)

    def tetragonal(self, directory="tmp-qe-relax-tetragonal", runopt="gen",  auto=0):
        """
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)

        shutil.copyfile(self.arts.xyz.file, os.path.join(directory, os.path.basename(self.arts.xyz.file)))
        #all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
        all_file = os.listdir()
        for element in self.arts.xyz.specie_labels:
            for item in all_file:
                if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                    shutil.copyfile(item, os.path.join(directory, item))
                    break
        self.arts.pseudo.dir = os.path.abspath(directory)
        self.control.pseudo_dir = os.path.abspath(directory)
        #
        os.chdir(directory)

        with open("relax.in.template", 'w') as fout:
            self.control.to_in(fout)
            self.system.to_in(fout)
            self.electrons.to_in(fout)
            self.ions.to_in(fout)

            coordtype = "angstrom"
            fout.write("ATOMIC_SPECIES\n")
            upf_all = [s for s in os.listdir(self.arts.pseudo.dir) if s.split(".")[-1] == "UPF"]
            for element in self.arts.xyz.specie_labels:
                for upf in upf_all:
                    if upf.split(".")[0] == element:
                        pseudo_file =upf
                        break
                fout.write("%s %f %s\n" % (element, base.element[element].mass, pseudo_file))
                pseudo_file = None
                # after pseudo_file used, set it to None to avoid it will be used in the next element
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

        # gen local bash script
        with open("relax-tetragonal.slurm", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
            fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
            fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
            fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
            fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
            fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
            #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

            a = self.arts.xyz.cell[0][0]
            c = self.arts.xyz.cell[2][2]

            fout.write("v11=%f\n" % self.arts.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.arts.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.arts.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.arts.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.arts.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.arts.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.arts.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.arts.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.arts.xyz.cell[2][2])

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("  for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("  do\n")
                    fout.write("    cp relax.in.template relax-${a}-${c}.in\n")
                    fout.write("    cat >> relax-${a}-${c}.in <<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("    yhrun $PMF_PWX < relax-${a}-${c}.in > relax-${a}-${c}.out\n")
                    fout.write("  done\n")
                else:
                    # only optimize a
                    fout.write("    cp relax.in.template relax-${a}.in\n")
                    fout.write("    cat >> relax-${a}.in <<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("    yhrun $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  cp relax.in.template relax-${c}.in\n")
                    fout.write("  cat >> relax-${c}.in<<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("0.000000 ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  yhrun $PMF_PWX < relax-${c}.in > relax-${c}.out\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass


        # gen pbs script
        with open("relax-tetragonal.pbs", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % self.run_params["jobname"])
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
            if "queue" in self.run_params and self.run_params["queue"] != None:
                fout.write("#PBS -q %s\n" %self.run_params["queue"])            
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

            a = self.arts.xyz.cell[0][0]
            c = self.arts.xyz.cell[2][2]

            fout.write("v11=%f\n" % self.arts.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.arts.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.arts.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.arts.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.arts.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.arts.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.arts.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.arts.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.arts.xyz.cell[2][2])

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("  for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("  do\n")
                    fout.write("    cp relax.in.template relax-${a}-${c}.in\n")
                    fout.write("    cat >> relax-${a}-${c}.in <<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("    mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < relax-${a}-${c}.in > relax-${a}-${c}.out\n")
                    fout.write("  done\n")
                else:
                    # only optimize a
                    fout.write("    cp relax.in.template relax-${a}.in\n")
                    fout.write("    cat >> relax-${a}.in <<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("    mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < relax-${a}.in > relax-${a}.out\n")
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  cp relax.in.template relax-${c}.in\n")
                    fout.write("  cat >> relax-${c}.in<<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("0.000000 ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < relax-${c}.in > relax-${c}.out\n")
                    fout.write("done\n")
                else:
                    # neither a or c is optimized
                    pass

        # gen local bash script
        with open("relax-tetragonal.bash", 'w') as fout:
            fout.write("#!/bin/bash\n")

            a = self.arts.xyz.cell[0][0]
            c = self.arts.xyz.cell[2][2]

            fout.write("v11=%f\n" % self.arts.xyz.cell[0][0])
            fout.write("v12=%f\n" % self.arts.xyz.cell[0][1])
            fout.write("v13=%f\n" % self.arts.xyz.cell[0][2])
            fout.write("v21=%f\n" % self.arts.xyz.cell[1][0])
            fout.write("v22=%f\n" % self.arts.xyz.cell[1][1])
            fout.write("v23=%f\n" % self.arts.xyz.cell[1][2])
            fout.write("v31=%f\n" % self.arts.xyz.cell[2][0])
            fout.write("v32=%f\n" % self.arts.xyz.cell[2][1])
            fout.write("v33=%f\n" % self.arts.xyz.cell[2][2])

            if na >= 2:
                # a is optimized
                fout.write("for a in `seq -w %f %f %f`\n" % (a-na/2*stepa, stepa, a+na/2*stepa))
                fout.write("do\n")
                if nc >= 2:
                    # optimize both a and c
                    fout.write("  for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("  do\n")
                    fout.write("    cp relax.in.template relax-${a}-${c}.in\n")
                    fout.write("    cat >> relax-${a}-${c}.in <<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("    %s pw.x < relax-${a}-${c}.in | tee relax-${a}-${c}.out\n" % self.run_params["mpi"])
                    fout.write("  done\n")
                else:
                    # only optimize a
                    fout.write("    cp relax.in.template relax-${a}.in\n")
                    fout.write("    cat >> relax-${a}.in <<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${a} 0.000000 0.000000\n")
                    fout.write("0.000000 ${a} 0.000000\n")
                    fout.write("0.000000 0.000000 ${v33}\n")
                    fout.write("EOF\n")
                    fout.write("    %s $PMF_PWX < relax-${a}.in | tee relax-${a}.out\n" % self.run_params["mpi"])
                fout.write("done\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only optimize c
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  cp relax.in.template relax-${c}.in\n")
                    fout.write("  cat >> relax-${c}.in<<EOF\n")
                    fout.write("\n")
                    fout.write("CELL_PARAMETERS angstrom\n")
                    fout.write("${v11} 0.000000 0.000000\n")
                    fout.write("0.000000 ${v22} 0.000000\n")
                    fout.write("0.000000 0.000000 ${c}\n")
                    fout.write("EOF\n")
                    fout.write("  %s $PMF_PWX < relax-${c}.in | tee relax-${c}.out\n" % self.run_params["mpi"])
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
                    fout.write("  for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
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
                    fout.write("set title Energy Latconst\n")
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
                    fout.write("set title Energy Latconst\n")
                    fout.write("set xlabel 'latconst(a)'\n")
                    fout.write("set ylabel 'Energy'\n")
                    fout.write("plot 'energy-latconst.data' w l\n")
                    fout.write("EOF\n")
            else:
                # a is not optimized
                if nc >= 2:
                    # only c is optimized
                    fout.write("for c in `seq -w %f %f %f`\n" % (c-nc/2*stepc, stepc, c+nc/2*stepc))
                    fout.write("do\n")
                    fout.write("  energy=`cat ../relax-${c}.out | grep '!    total energy' | tail -1`\n")
                    fout.write("  cat >> energy-latconst.data <<EOF\n")
                    fout.write("${c} ${energy:32:-2}\n")
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
                    # neither a nor c is optimized
                    pass
        #os.system("cd post-processing; bash get_energy.sh; cd ../")
        os.chdir("../")
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("bash relax-tetragonal.sh")
            os.chdir("../")
        server_handle(auto=args.auto, directory=directory, jobfilebase="relax-tetragonal", server=args.server)
