#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.qe.pwscf import pwscf


class static_run(pwscf):
    """
    About:
        static_run implements the control over scf, nscf and
        calculations based on them, like dos, pdos, bands, epsilon
    Status:
        currently implemented calculation including:
            scf, nscf, dos, bands, projwfc(pdos),
            epsilon, turbo_davidson, turbo_lanczos,
            phx_qmesh, q2r, matdyn, plotband, phx_gamma,
            dynmat, ir_raman, elf, fermi_surface,
            difference_charge_density, ellectron_density,

            converge test:
                ecutwfc, ecutrho, kpoints, degauss
    #
    """
    def __init__(self):
        super().__init__()
        
        self.control.basic_setting("scf") 

    def scf(self, directory="tmp-qe-static", inpname="static-scf.in", output="static-scf.out", mpi="", runopt="gen",
            jobname="pwscf-scf", nodes=1, ppn=32):
        """
        directory: a place for all the generated files

        parameters:
            directory: the overall static calculation directory

            runopt: determine whether the calculation is executed.
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
            shutil.copyfile(self.arts.xyz.file, os.path.join(directory, self.arts.xyz.file))
            all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
            for element in self.arts.xyz.specie_labels:
                for upf in all_upfs:
                    if upf.split(".")[0] == element:
                        shutil.copyfile(upf, os.path.join(directory, upf))
                        break
            # 

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)

            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="pw.x")
            # gen pbs scripts
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="pw.x", jobname=jobname, nodes=nodes, ppn=ppn)

        if runopt == 'genrun' or runopt == 'run':
            os.chdir(directory)
            os.system("%s pw.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def nscf(self, directory="tmp-qe-static", inpname="static-nscf.in", output="static-nscf.out", mpi="", runopt='gen',
            jobname="pscf-nscf", nodes=1, ppn=32):
        """
        parameters:
            directory: the overall static calculation directory

            runopt: determine whether the calculation is executed.
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

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)            
            
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="pw.x")
            # gen pbs scripts
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="pw.x", jobname=jobname, nodes=nodes, ppn=ppn)

        if runopt == 'genrun' or runopt == 'run':
            os.chdir(directory)
            os.system("%s pw.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")
    

    def converge_ecutwfc(self, emin, emax, step, directory="tmp-qe-ecutwfc", mpi="", runopt="gen",
            jobname="converge-ecutwfc", nodes=1, ppn=32):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.UPF %s/" % directory)
            os.system("cp %s %s/" % (self.arts.xyz.file, directory))
    
            os.chdir(directory)
            n_test = int((emax - emin) / step)
            for i in range(n_test + 1):
                ecut_wfc = int(emin + i * step)
                ecut_rho = ecut_wfc * 4 # using default value for ecut_rho: 4 * ecutwfc
                inp_name = "ecutwfc-%d.in" % ecut_wfc
                self.control.params['outdir'] = './tmp-' + str(ecut_wfc)
                self.system.params['ecutwfc'] = ecut_wfc
                with open(inp_name, 'w') as fout:
                    self.control.to_in(fout)
                    self.system.to_in(fout)
                    self.electrons.to_in(fout)
                    self.arts.to_in(fout)
            # gen yhbatch running script
            with open("converge-ecutwfc.sub", 'w') as fout:
                fout.write("#!/bin/bash\n")
                for i in range(n_test + 1):
                    ecut_wfc = int(emin + i * step)
                    inp_name = "ecutwfc-%d.in" % ecut_wfc
                    out_f_name = "ecutwfc-%d.out" % ecut_wfc
                    fout.write("yhrun -N 1 -n 24 pw.x < %s > %s\n" % (inp_name, out_f_name))
            # gen pbs running script
            with open("converge-ecutwfc.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % jobname)
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    ecut_wfc = int(emin + i * step)
                    inp_name = "ecutwfc-%d.in" % ecut_wfc
                    out_f_name = "ecutwfc-%d.out" % ecut_wfc
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE pw.x < %s > %s\n" % (inp_name, out_f_name))
            os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            for i in range(n_test + 1):
                ecut_wfc = int(emin + i * step)
                inp_name = "ecutwfc-%d.in" % ecut_wfc
                out_f_name = "ecutwfc-%d.out" % ecut_wfc
                os.system("%s pw.x < %s | tee %s" % (mpi, inp_name, out_f_name))
            os.chdir("../")

        
    def converge_ecutrho(self, emin, emax, step, ecutwfc, directory="tmp-qe-ecutrho", mpi="", runopt="gen",
            jobname="converge-ecutrho", nodes=1, ppn=32):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.UPF %s/" % directory)


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
            with open("converge-ecutrho.sub", 'w') as fout:
                fout.write("#!/bin/bash\n")
                for i in range(n_test + 1):
                    ecut_rho = int(emin + i * step)
                    inp_name = "ecutrho-%d.in" % ecut_rho
                    out_f_name = "ecutrho-%d.out" % ecut_rho
                    fout.write("yhrun -N 1 -n 24 pw.x < %s > %s\n" % (inp_name, out_f_name))
            # gen pbs running script
            with open("converge-ecutrho.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % jobname)
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    ecut_rho = int(emin + i * step)
                    inp_name = "ecutrho-%d.in" % ecut_rho
                    out_f_name = "ecutrho-%d.out" % ecut_rho
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE pw.x < %s > %s\n" % (inp_name, out_f_name))
            os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            for i in range(n_test + 1):
                ecut_rho = int(emin + i * step)
                inp_name = "ecutrho-%d.in" % ecut_rho
                out_f_name = "ecutrho-%d.out" % ecut_rho
                os.system("%s pw.x < %s | tee %s" % (mpi, inp_name, out_f_name))
            os.chdir("../")
    #
    def converge_kpoints(self, nk_min, nk_max, step=1, directory="tmp-qe-kpoints", mpi="", runopt="gen",
            jobname="converge-kpoints", nodes=1, ppn=32):
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
            os.system("cp *.UPF %s/" % directory)
	    
            os.chdir(directory)	
            n_test = int((nk_max - nk_min) / step)
            for i in range(n_test + 1):
                nk = nk_min + i * step # nk1 = nk2 = nk3 = nk
                inp_name = "kpoints-%d.in" % nk
                self.control.params['outdir'] = './tmp-' + str(nk)
                self.arts.set_kpoints([nk, nk, nk, 0, 0, 0])
                with open(inp_name, 'w') as fout:
                    self.control.to_in(fout)
                    self.system.to_in(fout)
                    self.electrons.to_in(fout)
                    self.arts.to_in(fout)
 
            # gen yhbatch running script
            with open("converge-kpoints.sub", 'w') as fout:
                fout.write("#!/bin/bash\n")
                for i in range(n_test + 1):
                    nk = nk_min + i * step # nk1 = nk2 = nk3 = nk
                    inp_name = "kpoints-%d.in" % nk
                    out_f_name = "kpoints-%d.out" % nk
                    fout.write("yhrun -N 1 -n 24 pw.x < %s > %s\n" % (inp_name, out_f_name))                   
            # gen pbs running script
            with open("converge-kpoints.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % jobname)
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    nk = nk_min + i * step # nk1 = nk2 = nk3 = nk
                    inp_name = "kpoints-%d.in" % nk
                    out_f_name = "kpoints-%d.out" % nk
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE pw.x < %s > %s\n" % (inp_name, out_f_name))
            os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            for i in range(n_test + 1):
                nk = nk_min + i * step
                inp_name = "kpoints-%d.in" % nk
                out_f_name = "kpoints-%d.out" % nk
                os.system("%s pw.x < %s | tee %s" % (mpi, inp_name, out_f_name))
            os.chdir("../")  

    def converge_degauss(self,degauss_min, degauss_max, step=0.01, directory="tmp-qe-degauss", mpi="",
            jobname="converge-degauss", nodes=1, ppn=32):
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
            os.system("cp *.UPF %s/" % directory)
	   
            os.chdir(directory)	
            n_test = int((degauss_max - degauss_min) / step)
            for i in range(n_test + 1):
                degauss = degauss_min + i * step
                inp_name = "degauss-%f.in" % degauss
                self.control.params['outdir'] = './tmp-' + str(degauss) 
                #self.arts.set_kpoints([nk, nk, nk, 0, 0, 0]) # use the previously convered kpoints(automatic)
                self.system.params['degauss'] = degauss
                with open(inp_name, 'w') as fout:
                    self.control.to_in(fout)
                    self.system.to_in(fout)
                    self.electrons.to_in(fout)
                    self.arts.to_in(fout)

            # gen yhbatch running script
            with open("converge-degauss.sub", 'w') as fout:
                fout.write("#!/bin/bash\n")
                for i in range(n_test + 1):
                    degauss = degauss_min + i * step
                    inp_name = "degauss-%f.in" % degauss
                    out_f_name = "degauss-%f.out" % degauss
                    fout.write("yhrun -N 1 -n 24 pw.x < %s > %s\n" % (inp_name, out_f_name))                   
            # gen pbs running script
            with open("converge-degauss.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % jobname)
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    degauss = degauss_min + i * step
                    inp_name = "degauss-%f.in" % degauss
                    out_f_name = "degauss-%f.out" % degauss
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE pw.x < %s > %s\n" % (inp_name, out_f_name))
            os.chdir("../")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            for i in range(n_test + 1):
                degauss = degauss_min + i * step
                inp_name = "degauss-%f.in" % degauss
                out_f_name = "degauss-%f.out" % degauss
                os.system("%s pw.x < %s | tee %s" % (mpi, inp_name, out_f_name))
            os.chdir("../")  

    
    def dos(self, directory="tmp-qe-static", inpname="static-dos.in", output="static-dos.out", mpi="",
            fildos="dosx.dos", bz_sum='smearing', ngauss='default', degauss='default', emin='default', emax='default',
            deltae='default', runopt="gen",
            jobname="dos", nodes=1, ppn=32):
        """
        Reference:
            http://www.quantum-espresso.org/Doc/INPUT_DOS.html
        
        bz_sum:
            'smearing' :
            'tetrahedra' :
            'tetrahedra_lin' :
            'tetrahedra_opt' :
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
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="dos.x")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="dos.x", jobname=jobname, nodes=nodes, ppn=ppn)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s dos.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def set_bands(self, bands_input={}):
        self.bands_input = {
                "prefix": self.control.params["prefix"],
                "outdir": self.control.params["outdir"],
                "filband": "bands.dat",
                "lsym": ".true."
                }
        for item in bands_input:
            self.bands_input[item] = bands_input[item]

    def bands(self, directory="tmp-qe-static", inpname1="static-bands.in", output1="static-bands.out",
            inpname2="bands.in", output2="bands.out", mpi="", runopt="gen", 
            jobname="band-structure", nodes=1, ppn=32):
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
            with open(os.path.join(directory, "band-structure.sub"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % ("pw.x", inpname1, output1))
                fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % ("bands.x", inpname2, output2))
            # gen pbs script
            with open(os.path.join(directory, "band-structure.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % jobname)
                fout.write("#PBS -l nodes=%d;ppn=%d\n" % (nodes, ppn))
                fout.write("mpirun -np %d -machinefile $PBS_NODEFILE %s < %s > %s\n" % (nodes*ppn, "pw.x", inpname1, output1))
                fout.write("mpirun -np %d -machinefile $PBS_NODEFILE %s < %s > %s\n" % (nodes*ppn, "bands.x", inpname2, output2))

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s pw.x < %s | tee %s" % (mpi, inpname1, output1))
            os.system("%s bands.x < %s | tee %s" % (mpi, inpname2, output2))
            os.chdir("../")

    def set_projwfc(self, projwfc_input={}):
        """
        """
        self.projwfc_input = {
                "prefix": self.control.params["prefix"],
                "outdir": self.control.params["outdir"],
                "filpdos": "projwfc",
                "ngauss": "default",
                "degauss": "default",
                "emin": "default",
                "emax": "default",
                "deltae": "default",
                }
        
        for item in projwfc_input:
            self.projwfc_input[item] = projwfc_input[item]


    def projwfc(self, directory="tmp-qe-static", inpname="static-projwfc.in", output="static-projwfc.out", mpi="", runopt="gen",
            jobname="projwfc-pdos", nodes=1, ppn=32):
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
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="projwfc.x")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="projwfc.x", jobname=jobname, nodes=nodes, ppn=ppn)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s projwfc.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")
    
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
            mpi="", runopt="gen", jobname="moledularpdos", nodes=1, ppn=32):
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
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="molecularpdos.x")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="molecularpdos.x", jobname=jobname, nodes=nodes, ppn=ppn)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s molecularpdos.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def epsilon(self, directory="tmp-qe-static", inpname="epsilon.in", output="epsilon.out", mpi="", runopt="gen",
            jobname="epsilon", nodes=1, ppn=32):
        """
        References:
            https://gitlab.com/QEF/material-for-ljubljana-qe-summer-school/blob/master/Day-3/handson-day3-TDDFPT.pdf

        epsilon.x:
            calculation of absorption spectra in IPA(Independent Particle 
            Approximation).
        Note:
            USPP rea not implemented now
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("epsilon.x calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("&INPUTPP\n")
                fout.write("calculation = 'eps'\n")
                fout.write("prefix = '%s'\n" % self.control.params["prefix"])
                fout.write("outdir = '%s'\n" % self.control.params["outdir"])
                fout.write("/\n")
                fout.write("\n")
                fout.write("&ENERGY_GRID\n")
                fout.write("smeartype = 'gaussian'\n") # type of smearing of the spectrum
                fout.write("intersmear = 0.1\n")       # the valus of smearing in eV
                fout.write("wmin = 0.0\n")             # minimum value of frequencies for a plot in eV
                fout.write("wmax = 15.0\n")            # maximum value of frequencies for a plot in eV
                fout.write("nw = 1000\n")              # number of points between wmin and wmax
                fout.write("/\n")
                fout.write("\n")

            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="epsilon.x")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="epsilon.x", jobname=jobname, nodes=nodes, ppn=ppn)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s epsilon.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def turbo_davidson(self, directory="tmp-qe-static", inpname1="turbo-davidson.in", output1="turbo-davidson.out",
            inpname2="turbo-spectrum-davidson.in", output2="turbo-spectrum-davidson.out", mpi="", runopt="gen",
            jobname="turbo-davidson", nodes=1, ppn=32):
        """
        References:
            https://gitlab.com/QEF/material-for-ljubljana-qe-summer-school/blob/master/Day-3/handson-day3-TDDFPT.pdf

        turbo_davidson.x:
            calculate adsorption spectra of molecules using time-dependent
            density functinal perturbation theory(TDDFPT).

            if if_dft_spectrum is set to .true. the result will be the same
            as calculated using epsilon.x, where there is no interaction.

            so set if_dft_spectrum to .false. to turn on the interaction.
            and you will get a shift of the peak compared to results calculated
            using epsilon.x or turbo_davidson.x(with if_dft_spectrum = .true.).

            when if_dft_spectrum is set to .true. turbo_davidson.x will produce
            a prefix-dft.eigen file, while a prefix.eigen file is produced if
            if_dft_spectrum is set to .false.

            we can also calculate absorption spectra using psudo-potential
            designed for B3LYP functional with turbo_davidson.x
            this way, we should set input_DFT = 'B3LYP' in scf calcualtion
            and set d0psi_rs = .true. in input file for turbo_davidson.x

        turbo_spectrum.x:
            post-processing calculation of the spectrum
        Note:
            turboTDDFT is not extended to metals, so we can only
            deal with insulators or semiconductors with turbo now.
            
            ltetra are not implemented now
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("turbo_davidson calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname1), 'w') as fout:
                fout.write("&lr_input\n")
                fout.write("prefix = '%s'\n" % self.control.params["prefix"])
                fout.write("outdir = '%s'\n" % self.control.params["outdir"])
                #fout.write("wfcdir = ''")
                fout.write("/\n")
                fout.write("\n")
                fout.write("&lr_dav\n")
                fout.write("if_dft_spectrum = .false.\n")
                fout.write("num_eign = 15\n")
                fout.write("num_init = 30\n")
                fout.write("num_basis_max = 90\n")
                fout.write("residue_conv_thr = 1.0E-6\n")
                fout.write("start = 0.0\n")
                fout.write("finish = 1.0\n")
                fout.write("step = 0.001\n")
                fout.write("broadening = 0.004\n")
                fout.write("reference = 0.3\n")
                fout.write("/\n")
                fout.write("\n")
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname1, output=output1, cmd="turbo_davidson.x")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname1, output=output1, cmd="turbo_davidson.x", jobname=jobname, nodes=nodes, ppn=ppn)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s turbo_davidson.x < %s | tee %s" % (mpi, inpname1, output1))
            os.chdir("../")
    
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname2), 'w') as fout:
                fout.write("&lr_input\n")
                fout.write("prefix = '%s'\n" % self.control.params["prefix"])
                fout.write("outdir = %s\n" % self.control.params["outdir"])
                fout.write("td = 'davidson'\n")
                fout.write("epsil = 0.004\n")
                fout.write("start = 0.0d0\n")
                fout.write("end = 1.0d0\n")
                fout.write("increment = 0.0001d0\n")
                fout.write("eign_file = 'pwscf.eigen'\n")
                fout.write("/\n")
                fout.write("\n")
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname2, output=output2, cmd="turbo_spectrum.x")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname2, output=output2, cmd="turbo_spectrum.x", jobname=jobname, nodes=nodes, ppn=ppn)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s turbo_spectrum.x < %s | tee %s" % (mpi, inpname2, output2))
            os.chdir("../")
    
    def turbo_lanczos(self, directory="tmp-qe-static", inpname1="turbo-lanczos.in", output1="turbo-lanczos.out",
            inpname2="turbo-spectrum-lanczos.in", output2="turbo-spectrum-lanczos.out", mpi="", runopt="gen",
            jobname="turbo_lanczos", nodes=1, ppn=32):
        """
        References:
            https://gitlab.com/QEF/material-for-ljubljana-qe-summer-school/blob/master/Day-3/handson-day3-TDDFPT.pdf

        turbo_lanczos.x:
            allows us to calculate absorption spectra of molecules using 
            time-dependent density functional perturbation theory (TDDFPT) 
            without computing empty states!
            
            turbo_lanczos.x allows us to obtain the absorption spectrum in a wide frequency
            range just by repeating a post-processing calculation using turbo_spectrum.x in a
            larger frequency range. This cannot be done with turbo_davidson.x

        turnbo_spectrum.x:
            post-processing calculation of the spectrum
            
        Note:
            turboTDDFT is not extended to metals
            ltetra are not implemented now
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("turbo_lanczos calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname1), 'w') as fout:
                fout.write("&lr_input\n")
                fout.write("prefix = '%s'\n" % self.control.params["prefix"])
                fout.write("outdir = '%s'\n" % self.control.params["outdir"])
                fout.write("restart_step = 100\n")
                fout.write("restart = .false.\n")
                fout.write("/\n")
                fout.write("\n")
                fout.write("&lr_control\n")
                fout.write("itermax = 1500\n")
                fout.write("ipol = 1\n")
                fout.write("/\n")
                fout.write("\n")
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname1, output=output1, cmd="turbo_lanczos.x")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname1, output=output1, cmd="turbo_lanczos.x", jobname=jobname, nodes=nodes, ppn=ppn)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s turbo_lanczos.x < %s | tee %s" % (mpi, inpname1, output1))
            os.chdir("../")
    
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname2), 'w') as fout:
                fout.write("&lr_input\n")
                fout.write("prefix = '%s'\n" % self.control.params["prefix"])
                fout.write("outdir = %s\n" % self.control.params["outdir"])
                fout.write("intermax0 = 1500\n") # Number of calculated Lanczos coefficient
                fout.write("intermax = 20000\n") # umber of extrapolated Lanczos coefficient
                fout.write("extrapolation = 'osc'\n") # Type of extrapolation (bi-constant)
                fout.write("epsil = 0.004\n") # The value of Lorenzian smearing in Ry
                fout.write("start = 0.0d0\n") # Minimum value of frequencies for a plot in Ry
                fout.write("end = 1.0d0\n") # Maximum value of frequencies for a plot in Ry
                fout.write("increment = 0.0001d0\n") # Frequency step in Ry
                fout.write("ipol = 1\n") # Polarization direction (same as in turbo_lanczos.x)
                fout.write("/\n")
                fout.write("\n")
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname2, output=output2, cmd="turbo_spectrum.x")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname2, output=output2, cmd="turbo_spectrum.x", jobname=jobname, nodes=nodes, ppn=ppn)
        
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s turbo_spectrum.x < %s | tee %s" % (mpi, inpname2, output2))
            os.chdir("../")


    def fermi_surface(self, directory="tmp-qe-static", inpname="fermi-surface.in", output="fermi-surface.out", mpi="", runopt="gen", jobname="fermi-surface", nodes=1, ppn=32):
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
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="fs.x")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="fs.x", jobname=jobname, nodes=nodes, ppn=ppn)
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s fs.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

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


    def pp(self, directory="tmp-qe-static", prefix="pp", mpi="", runopt="gen",
            jobname="pp.x-option", nodes=1, ppn=32):
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
            with open(os.path.join(directory, "pp.x.sub"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                for plot_num_i in self.inputpp["plot_num"]:
                    fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % ("pp.x", prefix+"-"+table[plot_num_i]+".in", prefix+"-"+table[plot_num_i]+".out"))   
            # gen pbs script
            with open(os.path.join(directory, "pp.x.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % jobname)
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
                for plot_num_i in self.inputpp["plot_num"]:
                    fout.write("mpirun -np %d -machinefile $PBS_NODEFILE %s < %s > %s\n" % (nodes*ppn, "pp.x", prefix+"-"+table[plot_num_i]+".in", prefix+"-"+table[plot_num_i]+".out"))   

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            for plot_num_i in self.inputpp["plot_num"]:
                os.system("%s pp.x < %s | tee %s" % (mpi, prefix+"-"+table[plot_num_i]+".in", prefix+"-"+table[plot_num_i]+".out"))
            os.chdir("../")


    def _pp_inputpp(self, fout, plot_num, filplot):
        """ 
        fout: a file stream for writing
        Note:
            plot_num -> selects what to save in filplot:
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
        fout: a file stream for writing
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

    def xspectra(self, directory="tmp-qe-static", inpname="xspectra.in", output="xspectra.out", mpi="", runopt="gen",
            jobname="xspectra", nodes=1, ppn=32):
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
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="xspectra.x")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="xspectra.x", jobname=jobname, nodes=nodes, ppn=ppn)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s xspectra.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")
    #

