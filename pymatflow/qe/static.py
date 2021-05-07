"""
Static calc
"""
import os
import re
import sys
import shutil


from pymatflow.remote.server import server_handle
from pymatflow.qe.pwscf import PwScf


class StaticRun(PwScf):
    """
    About:
        static_run implements the control over scf, nscf and
        calculations based on them, like dos, pdos, bands, epsilon
    Status:
        currently implemented calculation including:
        scf, nscf, dos, bands, projwfc(pdos), ir_raman, elf, fermi_surface,
        difference_charge_density, ellectron_density,
        converge test:
        ecutwfc, ecutrho, kpoints, degauss
    """
    def __init__(self):
        super().__init__()

        self.control.basic_setting("scf")
        self.arts.ifstatic = True

    def scf(self, directory="tmp-qe-static", inpname="static-scf.in", output="static-scf.out", runopt="gen", auto=0):
        """
        directory: a place for all the generated files

        :param directory: the overall static calculation directory

        :param runopt: determine whether the calculation is executed.
                there are three values: 'gen', 'genrun', 'run'
                'gen': only generate the input files
                'genrun': generate input files and run
                'run': run from the previously generated input files
        Note:
            only scf can generate the overall directory for static
            calculation(except the converge test for parameters like
            ecutwfc, kpoints, degauss)! other calculations is based
            on scf or nscf(which is based scf), so logically when
            doing these calculations there should already be the
            directory where scf calculation has been conducted.
        """
        self.control.calculation("scf")
        if runopt == 'gen' or runopt == 'genrun':
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
                self.arts.to_in(fout)

            # gen yhbatch script
            self.gen_llhpc(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX")
            # gen pbs scripts
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cdcloud script
            self.gen_cdcloud(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX")

        if runopt == 'genrun' or runopt == 'run':
            os.chdir(directory)
            os.system("%s $PMF_PWX < %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static-scf", server=self.run_params["server"])

    def nscf(self, directory="tmp-qe-static", inpname="static-nscf.in", output="static-nscf.out", runopt='gen', auto=0):
        """
        :param directory: the overall static calculation directory

        :param runopt: determine whether the calculation is executed.
                there are three values: 'gen', 'genrun', 'run'
                'gen': only generate the input files
                'genrun': generate input files and run
                'run': run from the previously generated input files
        """
        self.control.calculation("nscf")
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("non-scf calculation:\n")
            print("  directory of previous scf calculattion not found!\n")
            sys.exit(1)
        if runopt == 'gen' or runopt == 'genrun':

            self.arts.pseudo.dir = os.path.abspath(directory)
            self.control.set_params({"pseudo_dir": os.path.abspath(directory)})

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)

            # gen yhbatch script
            self.gen_llhpc(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX")
            # gen pbs scripts
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cdcloud script
            self.gen_cdcloud(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX")

        if runopt == 'genrun' or runopt == 'run':
            os.chdir(directory)
            os.system("%s $PMF_PWX < %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static-nscf", server=self.run_params["server"])

    def converge_ecutwfc(self, emin, emax, step, directory="tmp-qe-ecutwfc", runopt="gen", auto=0):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp %s %s/" % (self.arts.xyz.file, directory))
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

            os.chdir(directory)
            n_test = int((emax - emin) / step)
            for i in range(n_test + 1):
                ecut_wfc = int(emin + i * step)
                ecut_rho = ecut_wfc * 4 # using default value for ecut_rho: 4 * ecutwfc
                inp_name = "ecutwfc-%d.in" % ecut_wfc
                self.control.set_params({'outdir': './tmp-' + str(ecut_wfc)})
                self.system.set_params({'ecutwfc': ecut_wfc})
                with open(inp_name, 'w') as fout:
                    self.control.to_in(fout)
                    self.system.to_in(fout)
                    self.electrons.to_in(fout)
                    self.arts.to_in(fout)
            # gen yhbatch running script
            with open("converge-ecutwfc.slurm", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                for i in range(n_test + 1):
                    ecut_wfc = int(emin + i * step)
                    inp_name = "ecutwfc-%d.in" % ecut_wfc
                    out_f_name = "ecutwfc-%d.out" % ecut_wfc
                    fout.write("yhrun $PMF_PWX < %s > %s\n" % (inp_name, out_f_name))
            # gen pbs running script
            with open("converge-ecutwfc.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    ecut_wfc = int(emin + i * step)
                    inp_name = "ecutwfc-%d.in" % ecut_wfc
                    out_f_name = "ecutwfc-%d.out" % ecut_wfc
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < %s > %s\n" % (inp_name, out_f_name))
            os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            for i in range(n_test + 1):
                ecut_wfc = int(emin + i * step)
                inp_name = "ecutwfc-%d.in" % ecut_wfc
                out_f_name = "ecutwfc-%d.out" % ecut_wfc
                os.system("%s $PMF_PWX < %s | tee %s" % (self.run_params["mpi"], inp_name, out_f_name))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="converge-ecutwfc", server=self.run_params["server"])


    def converge_ecutrho(self, emin, emax, step, ecutwfc, directory="tmp-qe-ecutrho", runopt="gen", auto=0):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

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


            os.chdir(directory)
            n_test = int((emax - emin) / step)
            for i in range(n_test + 1):
                ecut_rho = int(emin + i * step)
                inp_name = "ecutrho-%d.in" % ecut_rho
                self.control.params['outdir'] = './tmp-' + str(ecut_rho)
                self.system.params['ecutwfc'] = ecutwfc
                self.system.params["ecutrho"] = ecut_rho
                with open(inp_name, 'w') as fout:
                    self.control.to_in(fout)
                    self.system.to_in(fout)
                    self.electrons.to_in(fout)
                    self.arts.to_in(fout)
            # gen yhbatch running script
            with open("converge-ecutrho.slurm", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                for i in range(n_test + 1):
                    ecut_rho = int(emin + i * step)
                    inp_name = "ecutrho-%d.in" % ecut_rho
                    out_f_name = "ecutrho-%d.out" % ecut_rho
                    fout.write("yhrun $PMF_PWX < %s > %s\n" % (inp_name, out_f_name))

            # gen pbs running script
            with open("converge-ecutrho.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    ecut_rho = int(emin + i * step)
                    inp_name = "ecutrho-%d.in" % ecut_rho
                    out_f_name = "ecutrho-%d.out" % ecut_rho
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < %s > %s\n" % (inp_name, out_f_name))
            os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            for i in range(n_test + 1):
                ecut_rho = int(emin + i * step)
                inp_name = "ecutrho-%d.in" % ecut_rho
                out_f_name = "ecutrho-%d.out" % ecut_rho
                os.system("%s $PMF_PWX < %s | tee %s" % (self.run_params["mpi"], inp_name, out_f_name))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="converge-ecutrho", server=self.run_params["server"])
    #
    def converge_kpoints(self, nk_min, nk_max, step=1, directory="tmp-qe-kpoints", runopt="gen", auto=0):
        """
        test the energy convergenc against k-points

        currently only support automatic schme of K_POINTS
        and only nk1 = nk2 = nk3 are supported

        Note:
            if you converge the ecutwfc previously, you should
            specify the converged ecutwfc through system in the
            parameters
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
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

            os.chdir(directory)
            n_test = int((nk_max - nk_min) / step)
            for i in range(n_test + 1):
                nk = nk_min + i * step # nk1 = nk2 = nk3 = nk
                inp_name = "kpoints-%d.in" % nk
                self.control.set_params({'outdir': './tmp-' + str(nk)})
                self.arts.set_kpoints([nk, nk, nk, 0, 0, 0])
                with open(inp_name, 'w') as fout:
                    self.control.to_in(fout)
                    self.system.to_in(fout)
                    self.electrons.to_in(fout)
                    self.arts.to_in(fout)

            # gen yhbatch running script
            with open("converge-kpoints.slurm", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                for i in range(n_test + 1):
                    nk = nk_min + i * step # nk1 = nk2 = nk3 = nk
                    inp_name = "kpoints-%d.in" % nk
                    out_f_name = "kpoints-%d.out" % nk
                    fout.write("yhrun $PMF_PWX < %s > %s\n" % (inp_name, out_f_name))
            # gen pbs running script
            with open("converge-kpoints.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    nk = nk_min + i * step # nk1 = nk2 = nk3 = nk
                    inp_name = "kpoints-%d.in" % nk
                    out_f_name = "kpoints-%d.out" % nk
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < %s > %s\n" % (inp_name, out_f_name))
            os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            for i in range(n_test + 1):
                nk = nk_min + i * step
                inp_name = "kpoints-%d.in" % nk
                out_f_name = "kpoints-%d.out" % nk
                os.system("%s $PMF_PWX < %s | tee %s" % (self.run_params["mpi"], inp_name, out_f_name))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="converge-kpoints", server=self.run_params["server"])

    def converge_degauss(self,degauss_min, degauss_max, step=0.01, directory="tmp-qe-degauss", runopt="gen", auto=0):
        """
        Convergence with respect to degauss/smearing

            smearing:
                (a) 'gauss'
                (b) 'marzari-vanderbilt'
                (c) 'methfessel-paxton'
            degauss:
                suggested values:
                0.06, 0.07, 0.08, 0.09, 0.10 (in Ry)
        Note:
            here we do the testing of degauss on energy.
            however quantities like the force on an atom
            may be more suited for this kind of testing.

            smearing is in fact part of the system setting
            how ever I set it a independent parameter in
            this function, to provide user the direct way
            to set the type of gauss smearing for testing.
            And of course we should not set smearing and
            occupations through system parameters.

            occpuations should always be set to smearing in
            testing degauss

            the user better set the previously converged
            ecutwfc throught system parameters
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
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

            os.chdir(directory)
            n_test = int((degauss_max - degauss_min) / step)
            for i in range(n_test + 1):
                degauss = degauss_min + i * step
                inp_name = "degauss-%f.in" % degauss
                self.control.set_params({'outdir': './tmp-%f' % degauss})
                #self.arts.set_kpoints([nk, nk, nk, 0, 0, 0]) # use the previously convered kpoints(automatic)
                self.system.set_params({'degauss': degauss})
                with open(inp_name, 'w') as fout:
                    self.control.to_in(fout)
                    self.system.to_in(fout)
                    self.electrons.to_in(fout)
                    self.arts.to_in(fout)

            # gen yhbatch running script
            with open("converge-degauss.slurm", 'w') as fout:
                fout.write("#!/bin/bash\n")
                for i in range(n_test + 1):
                    degauss = degauss_min + i * step
                    inp_name = "degauss-%f.in" % degauss
                    out_f_name = "degauss-%f.out" % degauss
                    fout.write("yhrun $PMF_PWX < %s > %s\n" % (inp_name, out_f_name))
            # gen pbs running script
            with open("converge-degauss.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    degauss = degauss_min + i * step
                    inp_name = "degauss-%f.in" % degauss
                    out_f_name = "degauss-%f.out" % degauss
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < %s > %s\n" % (inp_name, out_f_name))
            os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            for i in range(n_test + 1):
                degauss = degauss_min + i * step
                inp_name = "degauss-%f.in" % degauss
                out_f_name = "degauss-%f.out" % degauss
                os.system("%s $PMF_PWX < %s | tee %s" % (self.run_params["mpi"], inp_name, out_f_name))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="converge-degauss", server=self.run_params["server"])


    def dos(self, directory="tmp-qe-static", inpname="static-dos.in", output="static-dos.out",
            fildos="dosx.dos", bz_sum='smearing', ngauss='default', degauss='default', emin='default', emax='default',
            deltae='default', runopt="gen", auto=0):
        """
        Reference:
            http://www.quantum-espresso.org/Doc/INPUT_DOS.html

        :param bz_sum:
            'smearing' :
            'tetrahedra' :
            'tetrahedra_lin' :
            'tetrahedra_opt' :
        :param ngauss:
            'default': read from saved input for pw.x
                    0: Simple Gaussian (default)
                    1: Methfessel-Paxton of order 1
                   -1: Marzari-Vanderbilt "cold smearing"
                  -99: Fermi-Dirac function
        :param degauss:
            gaussian broadening, Ry (not eV!)
            'default':
            a floating number

        Note:
            the degauss in dos.x can significantly affect
            the  plotting of dos,
            but I don't know whether the degauss in scf
            and nscf also has such significant effect. if
            so, I might need provdie more ability to set
            appropriate degauss in scf and nscf running.
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("dos calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("&DOS\n")
                fout.write("prefix = '%s'\n" % self.control.params["prefix"])
                fout.write("outdir = '%s'\n" % self.control.params["outdir"])
                fout.write("fildos = '%s'\n" % fildos)
                #fout.write("bz_sum = '%s'\n" % bz_sum)
                if bz_sum == 'smearing':
                    if ngauss == 'default':
                        fout.write("! use ngauss read from input for pw.x store in xxx.save\n")
                    else:
                        fout.write("ngauss = %d\n" % ngauss)
                    if degauss == 'default':
                        fout.write("! use degauss read from input for pw.x stored in xxx.save\n")
                        fout.write("! or degauss = DeltaE, if DeltaE is specified\n")
                        fout.write("! we better set degauss and ngauss ourselves!\n")
                    else:
                        fout.write("degauss = %f\n" % degauss)
                if emin == 'default':
                    fout.write("!using default Emin: lower band value plus 3 times gauss smearing value\n")
                else:
                    fout.write("emin = %f\n" % emin)
                if emax == 'default':
                    fout.write("!using default Emax: upper band value minus 3 times gauss smearing value\n")
                else:
                    fout.write("emax = %f\n" % emax)
                if deltae == 'default':
                    fout.write("!using default DeltaE value\n")
                else:
                    fout.write("deltae = %f\n" % deltae)
                fout.write("/\n")
                fout.write("\n")

            # gen yhbatch script
            self.gen_llhpc(directory=directory, inpname=inpname, output=output, cmd="$PMF_DOSX")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_DOSX", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cdcloud script
            self.gen_cdcloud(directory=directory, inpname=inpname, output=output, cmd="$PMF_DOSX")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_DOSX < %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static-dos", server=self.run_params["server"])

    def set_bands(self, bands_input={}):
        self.bands_input = {
                "prefix": self.control.params["prefix"].as_val(t=str, dim=0),
                "outdir": self.control.params["outdir"].as_val(t=str, dim=0),
                "filband": "bands.dat",
                "lsym": ".true."
                }
        for item in bands_input:
            self.bands_input[item] = bands_input[item]

    def bands(self, directory="tmp-qe-static", inpname1="static-bands.in", output1="static-bands.out",
            inpname2="bands.in", output2="bands.out", runopt="gen", auto=0):
        """
        first check whether there is a previous scf running
        Note:
            the calculation of 'bands' is based on the previous scf or nscf running
            namely there must be the xxx.save/charge-density.dat for pw.x to read
            and do the bands calculation

        Warning:
            now we better use tpiba_b type kpoints setting!!! as only the postprocess of that
            type of band structure calculation is implemented now
        """
        self.control.calculation('bands')
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("bands calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)

        if runopt == "gen" or runopt == "genrun":

            self.arts.pseudo.dir = os.path.abspath(directory)
            self.control.pseudo_dir = os.path.abspath(directory)

            with open(os.path.join(directory, inpname1), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)
            #
            with open(os.path.join(directory, inpname2), 'w') as fout:
                fout.write("&BANDS\n")
                for item in self.bands_input:
                    if self.bands_input[item] is not None:
                        if type(self.bands_input[item]) == str and self.bands_input[item].lower() not in [".true.", ".false."]:
                            fout.write("%s = '%s'\n" % (item, self.bands_input[item]))
                        else:
                            fout.write("%s = %s\n" % (item, self.bands_input[item]))
                fout.write("/\n")
                fout.write("\n")

            # gen yhbatch script
            with open(os.path.join(directory, "band-structure.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                fout.write("yhrun %s < %s > %s\n" % ("$PMF_PWX", inpname1, output1))
                fout.write("yhrun %s < %s > %s\n" % ("$PMF_BANDSX", inpname2, output2))
            # gen pbs script
            with open(os.path.join(directory, "band-structure.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % ("$PMF_PWX", inpname1, output1))
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % ("$PMF_BANDSX", inpname2, output2))
            # gen cdcloud script
            with open(os.path.join(directory, "band-structure.slurm_cd"), 'w') as fout:
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
                fout.write("srun --mpi=pmix_v3 %s < %s > %s\n" % ("$PMF_PWX", inpname1, output1))
                fout.write("srun --mpi=pmix_v3 %s < %s > %s\n" % ("$PMF_BANDSX", inpname2, output2))

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_PWX < %s | tee %s" % (self.run_params["mpi"], inpname1, output1))
            os.system("%s $PMF_BANDSX < %s | tee %s" % (self.run_params["mpi"], inpname2, output2))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="band-structure.pbs", server=self.run_params["server"])

    def set_projwfc(self, projwfc_input={}):
        """
        """
        self.projwfc_input = {
                "prefix": self.control.params["prefix"].as_val(t=str, dim=0),
                "outdir": self.control.params["outdir"].as_val(t=str, dim=0),
                "filpdos": "projwfc",
                "ngauss": "default",
                "degauss": "default",
                "emin": "default",
                "emax": "default",
                "deltae": "default",
                }

        for item in projwfc_input:
            self.projwfc_input[item] = projwfc_input[item]


    def projwfc(self, directory="tmp-qe-static", inpname="static-projwfc.in", output="static-projwfc.out", runopt="gen", auto=0):
        """
        Reference:
            http://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html

        &projwfc can using projwfc.x to calculate Lowdin charges, spilling
        parameter, projected DOS

            ngauss:
                'default': read from saved input for pw.x
                    0: Simple Gaussian (default)
                    1: Methfessel-Paxton of order 1
                   -1: Marzari-Vanderbilt "cold smearing"
                  -99: Fermi-Dirac function
            degauss:
                gaussian broadening, Ry (not eV!)
                'default':
                a floating number

        Note:
            the degauss in projwfc.x can significantly affect
            the  plotting of dos,
            but I don't know whether the degauss in scf
            and nscf also has such significant effect. if
            so, I might need provdie more ability to set
            appropriate degauss in scf and nscf running.
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("projwfc calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)

        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("&PROJWFC\n")
                for item in self.projwfc_input:
                    if item in ["ngauss", "degauss", "emin", "emax", "deltae"]:
                        continue
                    if self.projwfc_input[item] is not None:
                        if type(self.projwfc_input[item]) == str:
                            fout.write("%s = '%s'\n" % (item, self.projwfc_input[item]))
                        else:
                            fout.write("%s = %s\n" % (item, self.projwfc_input[item]))
                if self.projwfc_input["ngauss"] == 'default':
                    fout.write("! use ngauss read from input for pw.x store in xxx.save\n")
                else:
                    fout.write("ngauss = %d\n" % self.projwfc_input["ngauss"])
                if self.projwfc_input["degauss"] == 'default':
                    fout.write("! use degauss read from input for pw.x stored in xxx.save\n")
                    fout.write("! or degauss = DeltaE, if DeltaE is specified\n")
                    fout.write("! we better set degauss and ngauss ourselves!\n")
                else:
                    fout.write("degauss = %f\n" % self.projwfc_input["degauss"])
                if self.projwfc_input["emin"] == 'default':
                    fout.write("!using default Emin: lower band value plus 3 times gauss smearing value\n")
                else:
                    fout.write("emin = %f\n" % self.projwfc_input["emin"])
                if self.projwfc_input["emax"] == 'default':
                    fout.write("!using default Emax: upper band value minus 3 times gauss smearing value\n")
                else:
                    fout.write("emax = %f\n" % self.projwfc_input["emax"])
                if self.projwfc_input["deltae"] == 'default':
                    fout.write("!using default DeltaE value\n")
                else:
                    fout.write("deltae = %f\n" % self.projwfc_input["deltae"])
                fout.write("/\n")
                fout.write("\n")

            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="$PMF_PROJWFCX")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_PROJWFCX", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cdcloud script
            self.gen_cdcloud(directory=directory, inpname=inpname, output=output, cmd="$PMF_PROJWFCX")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_PROJWFCX < %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static-projwfc.pbs", server=self.run_params["server"])

    def set_molecularpdos(self, inputmopdos={}):
        """
        Reference:
            http://www.quantum-espresso.org/Doc/INPUT_molecularpdos.html


            ngauss:
                    0: Simple Gaussian (default)
                    1: Methfessel-Paxton of order 1
                   -1: Marzari-Vanderbilt "cold smearing"
                  -99: Fermi-Dirac function
            degauss:
                gaussian broadening, Ry (not eV!)
                a floating number

        Note:
            I don't know why the run of molecularpdos.x in my computer is not stable
            with all the same condition, it sometimes run successfully, and when you
            execute again it might give 'STOP error reading file'. and when you again
            execute it, it might work!!! unbelievable
        """
        self.inputmopdos = {
                "fileout": "molecularpdos",
                "ngauss": 0,
                "degauss": 0.001,
                "emin": "default",
                "emax": "default",
                "deltae": "default",
                }
        for item in inputmopdos:
            if item in self.inputmopdos:
                self.inputmopdos[item] = inputmopdos[item]


    def molecularpdos(self, directory="tmp-qe-static", inpname="static-molecularpdos.in", output="static-molecularpdos.out",
            runopt="gen", auto=0):
        """
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("molecularpdos calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)

        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("&INPUTMOPDOS\n")
                fout.write("xmlfile_full = '%s'\n" % "./tmp/pwscf.save/atomic_proj.xml")
                fout.write("xmlfile_part = '%s'\n" % "./tmp/pwscf.save/atomic_proj.xml")
                fout.write("fileout = '%s'\n" % self.inputmopdos["fileout"])
                fout.write("ngauss = %d\n" % self.inputmopdos["ngauss"])
                fout.write("! default degauss is 0.0 which will calse float number erros\n")
                fout.write("! we better set degauss and ngauss ourselves!\n")
                fout.write("degauss = %f\n" % self.inputmopdos["degauss"])
                if self.inputmopdos["emin"] == 'default':
                    fout.write("!using default Emin: band extrema\n")
                else:
                    fout.write("emin = %f\n" % self.inputmopdos["emin"])
                if self.inputmopdos["emax"] == 'default':
                    fout.write("!using default Emax: band extrema\n")
                else:
                    fout.write("emax = %f\n" % self.inputmopdos["emax"])
                fout.write("!Note deltae is in unit of eV while other variables like degauss is Rydberg\n")
                if self.inputmopdos["deltae"] == 'default':
                    fout.write("!using default DeltaE value: 0.01 in unit of eV\n")
                else:
                    fout.write("deltae = %f\n" % self.inputmopdos["deltae"])
                fout.write("/\n")
                fout.write("\n")

            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="$PMF_MOLECULARPDOSX")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_MOLECULARPDOSX", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cdcloud script
            self.gen_cdcloud(directory=directory, inpname=inpname, output=output, cmd="$PMF_MOLECULARPDOSX")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_MOLECULARPDOSX < %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static-molecular-pdos", server=self.run_params["server"])


    def fermi_surface(self, directory="tmp-qe-static", inpname="fermi-surface.in", output="fermi-surface.out", runopt="gen", auto=0):
        """
        scf->nscf(with denser k points)->fs.x
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("fs.x calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("&fermi\n")
                fout.write("prefix = '%s'\n" % self.control.params["prefix"])
                fout.write("outdir = '%s'\n" % self.control.params["outdir"])
                fout.write("/\n")
                fout.write("\n")
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="$PMF_FSX")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_FSX", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cdcloud script
            self.gen_cdcloud(directory=directory, inpname=inpname, output=output, cmd="$PMF_FSX")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_FSX < %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="fermi-surface", server=self.run_params["server"])

    def set_pp(self, inputpp={}, plotpp={}):
        self.inputpp = {
                "plot_num": [0],
                }
        for item in inputpp:
            if item in self.inputpp:
                self.inputpp[item] = inputpp[item]

        self.plotpp = {
                "iflag": 3,
                "output_format": 5,
                }
        for item in plotpp:
            if item in self.plotpp:
                self.plotpp[item] = plotpp[item]


    def pp(self, directory="tmp-qe-static", prefix="pp", runopt="gen", auto=0):
        """
        Note:
            the 3D charge plot like electron localization function and charge density
            can be used to fabricate 2D plots using vesta software(Utilities/'2D Data Display').
            where you can set (hkl) and depth to plot.
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("pp.x calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            table = {
                    0: "electron-pseudo-charge-density",
                    1: "total-potential",
                    2: "local-ionic-potential",
                    3: "ldos",
                    4: "local-density-of-electronic-entropy",
                    5: "stm",
                    6: "spin-polar",
                    7: "molecular-orbitals",
                    8: "electron-local-function",
                    9: "charge-density-minus-superposition-of-atomic-densities",
                    10: "ILDOS",
                    11: "v_bare+v_H-potential",
                    12: "sawtooth-electric-field-potential",
                    13: "nocollinear-magnetization",
                    17: "all-electron-charge-density-paw-only",
                    18: "exchage-correlation-magnetic-field-noncollinear-case",
                    19: "reduced-density-gradient",
                    20: "product-of-charge-density-with-hessian",
                    21: "all-electron-density-paw-only",
                    }
            for plot_num_i in self.inputpp["plot_num"]:
                with open(os.path.join(directory, prefix+"-"+table[plot_num_i]+".in"), 'w') as fout:
                    self._pp_inputpp(fout, plot_num=plot_num_i, filplot=table[plot_num_i]+".dat")
                    self._pp_plot(fout, output_format=self.plotpp["output_format"], iflag=self.plotpp["iflag"], filepp=table[plot_num_i]+".dat")

            # gen yhbatch script
            with open(os.path.join(directory, "pp.x.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                for plot_num_i in self.inputpp["plot_num"]:
                    fout.write("yhrun %s < %s > %s\n" % ("$PMF_PPX", prefix+"-"+table[plot_num_i]+".in", prefix+"-"+table[plot_num_i]+".out"))
            # gen pbs script
            with open(os.path.join(directory, "pp.x.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))
                for plot_num_i in self.inputpp["plot_num"]:
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % ("$PMF_PPX", prefix+"-"+table[plot_num_i]+".in", prefix+"-"+table[plot_num_i]+".out"))
            # gen cdcloud script
            with open(os.path.join(directory, "pp.x.slurm_cd"), 'w') as fout:
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
                for plot_num_i in self.inputpp["plot_num"]:
                    fout.write("srun --mpi=pmix_v3 %s < %s > %s\n" % ("$PMF_PPX", prefix+"-"+table[plot_num_i]+".in", prefix+"-"+table[plot_num_i]+".out"))

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for plot_num_i in self.inputpp["plot_num"]:
                os.system("%s pp.x < %s | tee %s" % (self.run_params["mpi"], prefix+"-"+table[plot_num_i]+".in", prefix+"-"+table[plot_num_i]+".out"))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="pp.x", server=self.run_params["server"])


    def _pp_inputpp(self, fout, plot_num, filplot):
        """
        :param fout: a file stream for writing
        :param plot_num -> selects what to save in filplot:
             0  = electron (pseudo-)charge density
             1  = total potential V_bare + V_H + V_xc
             2  = local ionic potential V_bare
             3  = local density of states at specific energy or grid of energies
                  (number of states per volume, in bohr^3, per energy unit, in Ry)
             4  = local density of electronic entropy
             5  = STM images
                  Tersoff and Hamann, PRB 31, 805 (1985)
             8  = electron localization function (ELF)
             9  = charge density minus superposition of atomic densities
             13 = the noncollinear magnetization.

            About other value of plot_num, refere to the the input manual
            of pp.x:
                http://www.quantum-espresso.org/Doc/INPUT_PP.html

        """
        fout.write("&inputpp\n")
        fout.write("prefix = '%s'\n" % self.control.params["prefix"])
        fout.write("outdir = '%s'\n" % self.control.params["outdir"])
        fout.write("filplot = '%s'\n" % (filplot))
        fout.write("plot_num = %d\n" % plot_num)
        if plot_num == 0:
            fout.write("spin_component = %d\n" % 0)
        elif plot_num == 1:
            fout.write("spin_component = %d\n" % 0)
        elif plot_num == 3:
            pass
        elif plot_num == 5:
            pass
        elif plot_num == 7:
            fout.write("kpoint(1) = 1\n")
            fout.write("kpoint(2) = 2\n")
            fout.write("kband(1) = 1\n")
            fout.write('kband(2) = 2\n')
        elif plot_num == 10:
            pass
        elif plot_num == 17:
            pass
        fout.write("/\n")

    def _pp_plot(self, fout, filepp, iflag=3, output_format=5,
            e1=[2.0, 0.0, 0.0], e2=[0.0, 2.0, 0.0], e3=[0.0, 0.0, 2.0],
            x0=[0.0, 0.0, 0.0], nx=1000, ny=1000, nz=1000):
        """
        :param fout: a file stream for writing
        """
        #fout.write("&inputpp\n")
        #fout.write("/\n\n")
        fout.write("&plot\n")
        fout.write("nfile = 1\n")
        fout.write("filepp(1) = '%s'\n" % (filepp))
        fout.write("weight(1) = 1.0\n")
        fout.write("iflag = %d\n" % iflag)
        fout.write("output_format = %d\n" % output_format)
        if iflag == 0 or iflag == 1:
            fout.write("e1(1) = %f, e1(2) = %f, e1(3) = %f\n" % (e1[0], e1[1], e1[2]))
            fout.write("x0(1) = %f, x0(2) = %f, x0(3) = %f\n" % (x0[0], x0[1], x0[2]))
            fout.write("nx = %d\n" % nx)
        elif iflag == 2:
            fout.write("e1(1) = %f, e1(2) = %f, e1(3) = %f\n" % (e1[0], e1[1], e1[2]))
            fout.write("e2(1) = %f, e2(2) = %f, e2(3) = %f\n" % (e2[0], e2[1], e2[2]))
            fout.write("x0(1) = %f, x0(2) = %f, x0(3) = %f\n" % (x0[0], x0[1], x0[2]))
            fout.write("nx = %d, ny = %d\n" % (nx, ny))
        elif iflag == 3:
            fout.write("e1(1) = %f, e1(2) = %f, e1(3) = %f\n" % (e1[0], e1[1], e1[2]))
            fout.write("e2(1) = %f, e2(2) = %f, e2(3) = %f\n" % (e2[0], e2[1], e2[2]))
            fout.write("e3(1) = %f, e3(2) = %f, e3(3) = %f\n" % (e3[0], e3[1], e3[2]))
            fout.write("x0(1) = %f, x0(2) = %f, x0(3) = %f\n" % (x0[0], x0[1], x0[2]))
            fout.write("nx = %d, ny = %d, nz = %d\n" % (nx, ny, nz))
        elif iflag == 4:
            fout.write("radius = %f\n" % radius)
            fout.write("nx = %d, ny = %d\n" % (nx, ny))
        if output_format == 0:
            fout.write("fileout = '%s'\n" % (filepp.split(".")[0]+".1d.gp"))
        elif output_format == 2:
            fout.write("fileout = '%s'\n" % (filepp.split(".")[0]+".plotrho"))
        elif output_format == 3:
            fout.write("fileout = '%s'\n" % (filepp.split(".")[0]+".2d.xsf"))
        elif output_format == 5:
            fout.write("fileout = '%s'\n" % (filepp.split(".")[0]+".3d.xsf"))
        elif output_format == 6:
            fout.write("fileout = '%s'\n" % (filepp.split(".")[0]+".cube"))
        elif output_format == 7:
            fout.write("fileout = '%s'\n" % (filepp.split(".")[0]+"2d.gp"))
        fout.write("/\n")
        fout.write("\n")

    def xspectra(self, directory="tmp-qe-static", inpname="xspectra.in", output="xspectra.out", runopt="gen", auto=0):
        """
        Reference:
            http://www.quantum-espresso.org/Doc/INPUT_XSpectra.txt
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("xspectra.x calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("&input_xspectra\n")
                fout.write("/\n")
                fout.write("&plot\n")
                fout.write("/\n")
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="$PMF_XSPECTRAX")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_XSPECTRAX", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cdcloud script
            self.gen_cdcloud(directory=directory, inpname=inpname, output=output, cmd="$PMF_XSPECTRAX")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_XSPECTRAX < %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="xspectra", server=self.run_params["server"])
    #

    def run(self, directory="tmp-qe-static", runopt="gen", auto=0, kpath=None,
        kpoints_mp_scf=[1, 1, 1, 0, 0, 0], kpoints_mp_nscf=[3, 3, 3, 0, 0, 0]):
        """
        directory: a place for all the generated files

        :param directory: the overall static calculation directory

        :param runopt: determine whether the calculation is executed.
                there are three values: 'gen', 'genrun', 'run'
                'gen': only generate the input files
                'genrun': generate input files and run
                'run': run from the previously generated input files
        Note:
            scf, nscf, pdos, bands in a single run
        """
        if runopt == 'gen' or runopt == 'genrun':
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
            # check hybrid functional
            # in pw.x non-scf calc, hybrid functional is not allowed
            input_dft = self.system.params["input_dft"].as_val(t=str, dim=0) if self.system.params["input_dft"].as_val() is not None else None


            # 1) scf
            self.control.calculation("scf")
            self.set_kpoints(kpoints_option="automatic", kpoints_mp=kpoints_mp_scf)
            with open(os.path.join(directory, "static-scf.in"), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)

            # 2) nscf
            self.control.calculation("nscf")
            # hybrid functional calc is not allowed in non-scf calc
            self.set_kpoints(kpoints_option="automatic", kpoints_mp=kpoints_mp_nscf)
            if input_dft.lower() == "pbe0" or input_dft.lower() == "b3lyp" or input_dft.lower() == "hse":
                self.system.set_params({"input_dft": "pbe"})
            with open(os.path.join(directory, "static-nscf.in"), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)

            # 3) projwfc
            with open(os.path.join(directory, "static-projwfc.in"), 'w') as fout:
                fout.write("&PROJWFC\n")
                for item in self.projwfc_input:
                    if item in ["ngauss", "degauss", "emin", "emax", "deltae"]:
                        continue
                    if self.projwfc_input[item] is not None:
                        if type(self.projwfc_input[item]) == str:
                            fout.write("%s = '%s'\n" % (item, self.projwfc_input[item]))
                        else:
                            fout.write("%s = %s\n" % (item, self.projwfc_input[item]))
                if self.projwfc_input["ngauss"] == 'default':
                    fout.write("! use ngauss read from input for pw.x store in xxx.save\n")
                else:
                    fout.write("ngauss = %d\n" % self.projwfc_input["ngauss"])
                if self.projwfc_input["degauss"] == 'default':
                    fout.write("! use degauss read from input for pw.x stored in xxx.save\n")
                    fout.write("! or degauss = DeltaE, if DeltaE is specified\n")
                    fout.write("! we better set degauss and ngauss ourselves!\n")
                else:
                    fout.write("degauss = %f\n" % self.projwfc_input["degauss"])
                if self.projwfc_input["emin"] == 'default':
                    fout.write("!using default Emin: lower band value plus 3 times gauss smearing value\n")
                else:
                    fout.write("emin = %f\n" % self.projwfc_input["emin"])
                if self.projwfc_input["emax"] == 'default':
                    fout.write("!using default Emax: upper band value minus 3 times gauss smearing value\n")
                else:
                    fout.write("emax = %f\n" % self.projwfc_input["emax"])
                if self.projwfc_input["deltae"] == 'default':
                    fout.write("!using default DeltaE value\n")
                else:
                    fout.write("deltae = %f\n" % self.projwfc_input["deltae"])
                fout.write("/\n")
                fout.write("\n")

            # 4) band structure
            self.control.calculation('bands')
            self.set_kpoints(kpoints_option="crystal_b", crystal_b=kpath)
            if input_dft.lower() == "pbe0" or input_dft.lower() == "b3lyp" or input_dft.lower() == "hse":
                self.system.params["input_dft"] = "pbe"            
            with open(os.path.join(directory, "static-bands.in"), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)
            #
            with open(os.path.join(directory, "bands.in"), 'w') as fout:
                fout.write("&BANDS\n")
                for item in self.bands_input:
                    if self.bands_input[item] is not None:
                        if type(self.bands_input[item]) == str and self.bands_input[item].lower() not in [".true.", ".false."]:
                            fout.write("%s = '%s'\n" % (item, self.bands_input[item]))
                        else:
                            fout.write("%s = %s\n" % (item, self.bands_input[item]))
                fout.write("/\n")
                fout.write("\n")

            # 5) pp.x
            prefix="pp"
            table = {
                    0: "electron-pseudo-charge-density",
                    1: "total-potential",
                    2: "local-ionic-potential",
                    3: "ldos",
                    4: "local-density-of-electronic-entropy",
                    5: "stm",
                    6: "spin-polar",
                    7: "molecular-orbitals",
                    8: "electron-local-function",
                    9: "charge-density-minus-superposition-of-atomic-densities",
                    10: "ILDOS",
                    11: "v_bare+v_H-potential",
                    12: "sawtooth-electric-field-potential",
                    13: "nocollinear-magnetization",
                    17: "all-electron-charge-density-paw-only",
                    18: "exchage-correlation-magnetic-field-noncollinear-case",
                    19: "reduced-density-gradient",
                    20: "product-of-charge-density-with-hessian",
                    21: "all-electron-density-paw-only",
                    }
            for plot_num_i in self.inputpp["plot_num"]:
                with open(os.path.join(directory, prefix+"-"+table[plot_num_i]+".in"), 'w') as fout:
                    self._pp_inputpp(fout, plot_num=plot_num_i, filplot=table[plot_num_i]+".dat")
                    self._pp_plot(fout, output_format=self.plotpp["output_format"], iflag=self.plotpp["iflag"], filepp=table[plot_num_i]+".dat")

            # gen yhbatch script
            with open(os.path.join(directory, "static.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                fout.write("yhrun $PMF_PWX < static-scf.in > static-scf.out\n")
                fout.write("yhrun $PMF_PWX < static-nscf.in > static-nscf.out\n")
                fout.write("yhrun $PMF_PROJWFCX < static-projwfc.in > static-projwfc.out\n")
                fout.write("yhrun $PMF_PWX < static-bands.in > static-bands.out\n")
                fout.write("yhrun $PMF_BANDSX < bands.in > bands.out\n")
                for plot_num_i in self.inputpp["plot_num"]:
                    fout.write("yhrun %s < %s > %s\n" % ("$PMF_PPX", prefix+"-"+table[plot_num_i]+".in", prefix+"-"+table[plot_num_i]+".out"))

            # gen pbs script
            with open(os.path.join(directory, "static.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < static-scf.in > static-scf.out\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < static-nscf.in > static-nscf.out\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PROJWFCX < static-projwfc.in > static-projwfc.out\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_PWX < static-bands.in > static-bands.out\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_BANDSX < bands.in > bands.out\n")
                for plot_num_i in self.inputpp["plot_num"]:
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % ("$PMF_PPX", prefix+"-"+table[plot_num_i]+".in", prefix+"-"+table[plot_num_i]+".out"))

            # gen local bash script
            with open(os.path.join(directory, "static.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("%s $PMF_PWX < static-scf.in | tee static-scf.out\n" % self.run_params["mpi"])
                fout.write("%s $PMF_PWX < static-nscf.in | tee static-nscf.out\n" % self.run_params["mpi"])
                fout.write("%s $PMF_PROJWFCX < static-projwfc.in | tee static-projwfc.out\n" % self.run_params["mpi"])
                fout.write("%s $PMF_PWX < static-bands.in > static-bands.out\n" % self.run_params["mpi"])
                fout.write("%s $PMF_BANDSX < bands.in | tee bands.out\n" % self.run_params["mpi"])
                for plot_num_i in self.inputpp["plot_num"]:
                    fout.write("%s $PMF_PPX < %s | tee %s\n" % (self.run_params["mpi"], prefix+"-"+table[plot_num_i]+".in", prefix+"-"+table[plot_num_i]+".out"))

            # gen cdcloud script
            with open(os.path.join(directory, "static.slurm_cd"), 'w') as fout:
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
                fout.write("srun --mpi=pmix_v3 $PMF_PWX < static-scf.in > static-scf.out\n")
                fout.write("srun --mpi=pmix_v3 $PMF_PWX < static-nscf.in > static-nscf.out\n")
                fout.write("srun --mpi=pmix_v3 $PMF_PROJWFCX < static-projwfc.in > static-projwfc.out\n")
                fout.write("srun --mpi=pmix_v3 $PMF_PWX < static-bands.in > static-bands.out\n")
                fout.write("srun --mpi=pmix_v3 $PMF_BANDSX < bands.in > bands.out\n")
                for plot_num_i in self.inputpp["plot_num"]:
                    fout.write("srun --mpi=pmix_v3 %s < %s > %s\n" % ("$PMF_PPX", prefix+"-"+table[plot_num_i]+".in", prefix+"-"+table[plot_num_i]+".out"))

        if runopt == 'genrun' or runopt == 'run':
            os.chdir(directory)
            os.system("bash static.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static", server=self.run_params["server"])
