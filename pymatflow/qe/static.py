#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.qe.base.control import qe_control
from pymatflow.qe.base.system import qe_system
from pymatflow.qe.base.electrons import qe_electrons
from pymatflow.qe.base.arts import qe_arts


class static_run:
    """
    About:
        static_run implements the control over scf, nscf and
        calculations based on them, like dos, pdos, bands, epsilon
        Phonon(DFPT)
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
    def __init__(self, xyz_f):
        self.control = qe_control()
        self.system = qe_system()
        self.electrons = qe_electrons()
        self.arts = qe_arts(xyz_f)

        self.control.basic_setting("scf") 
        self.system.basic_setting(self.arts)
        self.electrons.basic_setting()
        self.arts.basic_setting(ifstatic=True)


    def scf(self, directory="tmp-qe-static", inpname="static-scf.in", output="static-scf.out", 
            mpi="", runopt="gen", control={}, system={}, electrons={}, kpoints_option="automatic", kpoints_mp=[2, 2, 2, 0, 0, 0]):
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
        if runopt == 'gen' or runopt == 'genrun':
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.UPF %s/" % directory)
            os.system("cp %s %s/" % (self.arts.xyz.file, directory))
            
            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.set_occupations() which uses
            # self.system.set_occupations() to set them, as self.system.set_params() 
            # is suppressed from setting occupations related parameters
            self.set_occupations(system)
            self.control.set_params(control)
            self.system.set_params(system)
            self.electrons.set_params(electrons)
            self.arts.set_kpoints(option=kpoints_option, kpoints_mp=kpoints_mp)
            self.control.calculation("scf")
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="pw.x")

        if runopt == 'genrun' or runopt == 'run':
            os.chdir(directory)
            os.system("%s pw.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def nscf(self, directory="tmp-qe-static", inpname="static-nscf.in", output="static-nscf.out",
            mpi="", runopt='gen', control={}, system={}, electrons={}, kpoints_mp=[4, 4, 4, 0, 0, 0]):
        """
        parameters:
            directory: the overall static calculation directory

            runopt: determine whether the calculation is executed.
                there are three values: 'gen', 'genrun', 'run'
                'gen': only generate the input files
                'genrun': generate input files and run
                'run': run from the previously generated input files
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("non-scf calculation:\n")
            print("  directory of previous scf calculattion not found!\n")
            sys.exit(1)
        if runopt == 'gen' or runopt == 'genrun':

            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.set_occupations() which uses
            # self.system.set_occupations() to set them, as self.system.set_params() 
            # is suppressed from setting occupations related parameters
            self.set_occupations(system)
            self.control.set_params(control)
            self.system.set_params(system)
            self.electrons.set_params(electrons)
            self.arts.set_kpoints(kpoints_mp)
            self.control.calculation("nscf")
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)            
            
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="pw.x")

        if runopt == 'genrun' or runopt == 'run':
            os.chdir(directory)
            os.system("%s pw.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")
    

    def set_occupations(self, system):
        """
            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.system.set_occupations() to 
            # set them, as self.system.set_params() is suppressed from setting
            # occupations related parameters
            # if occupations == None, use default smearing occupation. and 
            # if occupations == "tetrahedra" the value set for smearing and degauss is ignored.
            # if occupations == "smearing", the value of smearing and degauss
            # should be legacy, not None or other illegal values.
        """
        if "occupations" in system:
            if system["occupations"] == None: # user default setting of set_occupations()
                self.system.set_occupations()
            elif system["occupations"] == "tetrahedra":
                self.system.set_occupations(occupations="tetrahedra")
            elif system["occupations"] == "smearing":
                if "smearing" in system and "degauss" in system:
                    self.system.set_occupations(occupations="smearing", smearing=system["smearing"], degauss=system["degauss"])
                elif "smearing" in system:
                    self.system.set_occupations(occupations="smearing", smearing=system["smearing"])
                elif "degauss" in system:
                    self.system.set_occupations(occupations="smearing", degauss=system["degauss"])
                else:
                    self.system.set_occupations(occupations="smearing")
            else:
                pass

    def converge_ecutwfc(self, emin, emax, step, directory="tmp-qe-ecutwfc", 
            mpi="", runopt="gen", control={}, system={}, electrons={}, kpoints_option="automatic", kpoints_mp=[1, 1, 1, 0, 0, 0]):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.UPF %s/" % directory)
            os.system("cp %s %s/" % (self.arts.xyz.file, directory))
    
            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.set_occupations() which uses
            # self.system.set_occupations() to set them, as self.system.set_params() 
            # is suppressed from setting occupations related parameters
            self.set_occupations(system)
            self.control.set_params(control)
            self.system.set_params(system)
            self.electrons.set_params(electrons)
            self.arts.set_kpoints(option=kpoints_option, kpoints_mp=kpoints_mp)
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

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            for i in range(n_test + 1):
                ecut_wfc = int(emin + i * step)
                inp_name = "ecutwfc-%d.in" % ecut_wfc
                out_f_name = "ecutwfc-%d.out" % ecut_wfc
                os.system("%s pw.x < %s | tee %s" % (mpi, inp_name, out_f_name))

            # analyse the result
            for i in range(n_test + 1):
                ecut_wfc = int(emin + i * step)
                out_f_name = "ecutwfc-%d.out" % ecut_wfc
                os.system("cat %s | grep '!    total energy' >> energy-ecutwfc.data" % out_f_name)

            ecut_wfc_all = [emin + i * step for i in range(n_test + 1)]
            energy_all = []
            with open("energy-ecutwfc.data", 'r') as fin:
                for line in fin:
                    energy_all.append(float(line.split()[4]))

            plt.plot(ecut_wfc_all, energy_all, marker='o')
            plt.title("Ecutwfc Converge Test", fontweight='bold', color='red')
            plt.xlabel("Ecutwfc (Ry)")
            plt.ylabel("Energy (Ry)")
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-ecutwfc.png")
            plt.show()
            os.chdir("../")

        
    def converge_ecutrho(self, emin, emax, step, ecutwfc, directory="tmp-qe-ecutrho", 
            mpi="", runopt="gen", control={}, system={}, electrons={}, kpoints_option="automatic", kpoints_mp=[1, 1, 1, 0, 0, 0]):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.UPF %s/" % directory)

            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.set_occupations() which uses
            # self.system.set_occupations() to set them, as self.system.set_params() 
            # is suppressed from setting occupations related parameters
            self.set_occupations(system)
            self.control.set_params(control)
            self.system.set_params(system)
            self.electrons.set_params(electrons)
            self.arts.set_kpoints(option=kpoints_option, kpoints_mp=kpoints_mp)

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

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            for i in range(n_test + 1):
                ecut_rho = int(emin + i * step)
                inp_name = "ecutrho-%d.in" % ecut_rho
                out_f_name = "ecutrho-%d.out" % ecut_rho
                os.system("%s pw.x < %s | tee %s" % (mpi, inp_name, out_f_name))
            # analyse the result
            for i in range(n_test + 1):
                ecut_rho = int(emin + i * step)
                out_f_name = "ecutrho-%d.out" % ecut_rho
                os.system("cat %s | grep '!    total energy' >> energy-ecutrho.data" % out_f_name)

            ecut_rho_all = [emin + i * step for i in range(n_test + 1)]
            energy_all = []
            with open("energy-ecutrho.data", 'r') as fin:
                for line in fin:
                    energy_all.append(float(line.split()[4]))

            plt.plot(ecut_rho_all, energy_all, marker='o')
            plt.title("Ecutrho Converge Test", fontweight='bold', color='red')
            plt.xlabel("Ecutrho (Ry)")
            plt.ylabel("Energy (Ry)")
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-ecutrho.png")
            plt.show()
            os.chdir("../")
    #
    def converge_kpoints(self, nk_min, nk_max, step=1, directory="tmp-qe-kpoints", 
            mpi="", control={}, system={}, electrons={}, runopt="gen"):
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
	    
            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.set_occupations() which uses
            # self.system.set_occupations() to set them, as self.system.set_params() 
            # is suppressed from setting occupations related parameters
            self.set_occupations(system)
            self.control.set_params(control)
            self.system.set_params(system)
            self.electrons.set_params(electrons)
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

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            for i in range(n_test + 1):
                nk = nk_min + i * step
                inp_name = "kpoints-%d.in" % nk
                out_f_name = "kpoints-%d.out" % nk
                os.system("%s pw.x < %s | tee %s" % (mpi, inp_name, out_f_name))

            # analyse the result
            for i in range(n_test + 1):
                nk = nk_min + i * step
                out_f_name = "kpoints-%d.out" % nk
                os.system("cat %s | grep '!    total energy' >> energy-kpoints.data" % out_f_name)
            
            nk_all = [nk_min + i * step for i in range(n_test + 1)]
            energy_all = []
            with open("energy-kpoints.data", 'r') as fin:
                for line in fin:
                    energy_all.append(float(line.split()[4]))
            
            plt.plot(nk_all, energy_all, marker='o')
            plt.title("kpoints converge test", fontweight='bold', color='red')
            plt.xlabel("Kpoints")
            plt.ylabel("Energy (Ry)")
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-kpoints.png")
            plt.show()
            os.chdir("../")  

    def converge_degauss(self,degauss_min, degauss_max, step=0.01, directory="tmp-qe-degauss", mpi="",
            smearing='gauss', control={}, system={}, electrons={}, runopt="gen", kpoints_mp=[1, 1, 1, 0, 0, 0]):
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
	   
            self.control.set_params(control)
            self.system.set_params(system)
            self.electrons.set_params(electrons)
            self.arts.set_kpoints(kpoints_mp)
            self.system.params["occupations"] = "smearing"

            os.chdir(directory)	
            n_test = int((degauss_max - degauss_min) / step)
            for i in range(n_test + 1):
                degauss = degauss_min + i * step
                inp_name = "degauss-%f.in" % degauss
                self.control.params['outdir'] = './tmp-' + str(degauss) 
                #self.arts.set_kpoints([nk, nk, nk, 0, 0, 0]) # use the previously convered kpoints(automatic)
                self.system.params['smearing'] = smearing
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
        if runopt == "run" or runopt == "genrun":
            # run the simulation
            for i in range(n_test + 1):
                degauss = degauss_min + i * step
                inp_name = "degauss-%f.in" % degauss
                out_f_name = "degauss-%f.out" % degauss
                os.system("%s pw.x < %s | tee %s" % (mpi, inp_name, out_f_name))

            # analyse the result
            for i in range(n_test + 1):
                degauss = degauss_min + i * step
                out_f_name = "degauss-%f.out" % degauss
                os.system("cat %s | grep '!    total energy' >> energy-degauss.data" % out_f_name)

            degauss_all = [degauss_min + i * step for i in range(n_test + 1)]
            energy_all = []
            with open("energy-degauss.data", 'r') as fin:
                for line in fin:
                    energy_all.append(float(line.split()[4]))

            plt.plot(degauss_all, energy_all, marker='o')
            plt.title("degauss converge test", fontweight='bold', color='red')
            plt.xlabel("degauss")
            plt.ylabel("Energy (Ry)")
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-degauss.png")
            plt.show()
            os.chdir("../")  

    
    def dos(self, directory="tmp-qe-static", inpname="static-dos.in", output="static-dos.out", mpi="",
            fildos="dosx.dos", bz_sum='smearing', ngauss='default', degauss='default', emin='default', emax='default',
            deltae='default', runopt="gen"):
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

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s dos.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def bands(self, directory="tmp-qe-static", inpname1="static-bands.in", output1="static-bands.out",
            inpname2="bands.in", output2="bands.out", mpi="", kptopt="automatic",
            control={}, system={}, electrons={}, kpoints_mp=[4, 4, 4, 0, 0, 0], runopt="gen"
            ):
        """
        first check whether there is a previous scf running
        Note:
            the calculation of 'bands' is based on the previous scf or nscf running
            namely there must be the xxx.save/charge-density.dat for pw.x to read
            and do the bands calculation
        """
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("bands calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        
        if runopt == "gen" or runopt == "genrun":
            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.set_occupations() which uses
            # self.system.set_occupations() to set them, as self.system.set_params() 
            # is suppressed from setting occupations related parameters
            self.set_occupations(system)
            self.control.set_params(control)
            self.system.set_params(system)
            self.electrons.set_params(electrons)
            self.control.calculation('bands')
            # ===========
            # set kpoints
            # ===========
            if kptopt == "automatic":
                self.arts.set_kpoints(kpoints_mp)
            elif kptopt == "crystal_b":
                self.arts.set_kpoints(option="crystal_b")

            with open(os.path.join(directory, inpname1), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)

            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="pw.x")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s pw.x < %s | tee %s" % (mpi, inpname1, output1))
            os.chdir("../")

        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname2), 'w') as fout:
                fout.write("&BANDS\n")
                fout.write("prefix = '%s'\n" % self.control.params["prefix"])
                fout.write("outdir = '%s'\n" % ("./tmp"))
                fout.write("filband = '%s'\n" % ("bands.dat"))
                fout.write("lsym = %s\n" % (".true."))
                fout.write("/\n")
                fout.write("\n")
            
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="bands.x")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s bands.x < %s | tee %s" % (mpi, inpname2, output2))
            os.chdir("../")
        

    def projwfc(self, directory="tmp-qe-static", inpname="static-projwfc.in", output="static-projwfc.out", 
            mpi="", filpdos="projwfc", ngauss='default', degauss='default', emin='default', emax='default',
            deltae='default', runopt="gen"
            ):
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
                fout.write("prefix = '%s'\n" % self.control.params["prefix"])
                fout.write("outdir = '%s'\n" % self.control.params["outdir"])
                fout.write("filpdos = '%s'\n" % filpdos)
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
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="projwfc.x")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s projwfc.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def epsilon(self, directory="tmp-qe-static", inpname="epsilon.in", output="epsilon.out", mpi="", runopt="gen"):
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

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s epsilon.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def turbo_davidson(self, directory="tmp-qe-static", inpname1="turbo-davidson.in", output1="turbo-davidson.out",
            inpname2="turbo-spectrum-davidson.in", output2="turbo-spectrum-davidson.out", mpi="", runopt="gen"):
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

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s turbo_spectrum.x < %s | tee %s" % (mpi, inpname2, output2))
            os.chdir("../")
    
    def turbo_lanczos(self, directory="tmp-qe-static", inpname1="turbo-lanczos.in", output1="turbo-lanczos.out",
            inpname2="turbo-spectrum-lanczos.in", output2="turbo-spectrum-lanczos.out", mpi="", runopt="gen"):
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
        
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s turbo_spectrum.x < %s | tee %s" % (mpi, inpname2, output2))
            os.chdir("../")

    def phx_qmesh(self, directory="tmp-qe-static", inpname="phx-qmesh.in", output="phx-qmesh.out", 
            dynamat_file="phx-qmesh.dyn", mpi="", runopt="gen", qpoints=[2, 2, 2]):
        """
        Reference:
            https://gitlab.com/QEF/material-for-ljubljana-qe-summer-school/blob/master/Day-3/handson-day3-DFPT.pdf
            http://www.quantum-espresso.org/Doc/ph_user_guide/
            http://www.fisica.uniud.it/~giannozz/QE-Tutorial/handson_phon.html

        ph.x:
            performing phonon calculation based on scf using 
            DFPT theory. it is the executable of PHonon package
            if qe.
            
            parameter epsil = .true. will calculate and store the
            dielectric tensor and effective charges which is for
            the polar materials

            we can do phonon calculation only at \Gamma point and also
            at a q-grid to get a phonon dispersion graph
            
        Note:
            this function(phonon_phx_qmesh) do phonon calculation at 
            a q-grid
            another function(phonon_phx_gamma) do phonon calculation only
            at \Gamma point

            PHonon: linear-response calculations(phonons, dielectric properties)
                (1) phonon frequencies and eigenvectors at a generic wave vector
                (2) dielectric tensor, effective charges, IR cross sections
                (3) interatomic force constants in real space
                (4) electron-phonon interaction coefficients for metals
                (5) nonresonant Raman cross sections
                (6) third-order anharmonic phonon lifetimes cross sections

        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("ph.x with qmesh calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("&inputph\n")
                fout.write("tr2_ph = 1.0d-14\n")
                fout.write("ldisp = .true.\n") # option for the calculation on a grid
                fout.write("nq1 = %d\n" % qpoints[0]) # 4
                fout.write("nq2 = %d\n" % qpoints[1]) # 4
                fout.write("nq3 = %d\n" % qpoints[2]) # 4
                fout.write("prefix = '%s'\n" % self.control.params["prefix"])
                fout.write("outdir = '%s'\n" % self.control.params["outdir"])
                fout.write("fildyn = '%s'\n" % dynamat_file)
                fout.write("/\n")
                fout.write("\n")
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="ph.x")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s ph.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def q2r(self, directory="tmp-qe-static", inpname="q2r.in", output="q2r.out", 
            dynamat_file="phx-qmesh.dyn", ifc_file="ifc.fc", mpi="", runopt="gen", zasr='simple'):
        """
        q2r.x:
            calculation of Interatomic Force Constants(IFC) from 
            Dynamical matrices from the phonon calculation
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("q2r calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("&input\n")
                fout.write("fildyn = '%s'\n" % dynamat_file) # Dynamical matrices from the phonon calculation
                fout.write("zasr = '%s'\n" % zasr) # A way to impose the acoustic sum rule
                fout.write("flfrc = '%s'\n" % ifc_file) # Output file of the interatomic force constants
                fout.write("/\n")
                fout.write("\n")
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="q2r.x")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s q2r.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def matdyn(self, directory="tmp-qe-static", inpname="matdyn.in", output="matdyn.out", 
            ifc_file="ifc.fc", mpi="", runopt="gen", asr='simple', 
            nqpoints=2, qpoints = [[0.0, 0.0, 0.0, 0.0], [0.012658, 0.0, 0.0, 0.012658]]):
        """
        matdyn.x
            Calculate phonons at generic q points using IFC
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("matdyn.x calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("&input\n")
                fout.write("asr = '%s'\n" % asr) # Acoustic sum rule
                fout.write("flfrc = '%s'\n" % ifc_file) # File with IFC's
                fout.write("flfrq = 'frequencies.freq'\n") # Output file with the frequencies
                fout.write("/\n")
                fout.write("%d\n" % nqpoints) # Number of q points
                for i in range(nqpoints):
                    fout.write("%f %f %f %f\n" % (qpoints[i][0], qpoints[i][1], qpoints[i][2], qpoints[i][3]))
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="matdyn.x")
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s matdyn.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def plotband(self, directory="tmp-qe-static", inpname="plotband.in", output="plotband.out", frequencies_file="frequencies.freq", mpi="", runopt="gen", freq_min=0, freq_max=600):
        """
        plotband.x
            Plot the phonon dispersion
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("plotband calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("%s\n" % frequencies_file) # Input file with the frequencies at various q
                fout.write("%f %f\n" % (freq_min, freq_max)) # Range of frequencies for a visualization
                fout.write("freq.plot\n") # Output file with frequencies which will be used for plot
                fout.write("freq.ps\n") # Plot of the dispersion
                fout.write("0.0\n") # Fermi level (needed only for band structure plot)
                fout.write("100.0 0.0\n") # Freq. step and reference freq. on the plot freq.ps
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="plotband.x")
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s plotband.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")


    def phx_gamma(self, directory="tmp-qe-static", inpname="phx-gamma.in", output="phx-gamma.out", 
            dynamat_file="phx-gamma.dyn", mpi="", runopt="gen"):
        """
        do phonon calculation only at \Gamma point
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("ph.x with gamma calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("&inputph\n")
                fout.write("tr2_ph = 1.0d-14\n")
                fout.write("prefix = '%s'\n" % self.control.params["prefix"])
                fout.write("outdir = '%s'\n" % self.control.params["outdir"])
                fout.write("fildyn = '%s'\n" % dynamat_file)
                fout.write("/\n")
                fout.write("0.0 0.0 0.0\n")
                fout.write("\n")
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="ph.x")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s ph.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def dynmat(self, directory="tmp-qe-static", inpname="dynmat.in", output="dynmat.out", 
            dynamat_file="phx-qmesh.dyn", mpi="", runopt="gen"):
        """
        imposing acoustic sum rule (ASR)
        extract the phonon information from ph.x output using dynmat.x(
        which can also be used to get IR and Raman.
        )
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("dynmat.x calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("&input\n")
                fout.write("fildyn = '%s'\n" % dynamat_file) # File containing the dynamical matrix
                fout.write("asr = 'simple'\n")
                fout.write("/\n")
                fout.write("\n")
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="dynmat.x")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s dynmat.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def ir_raman(self, directory="tmp-qe-static", mpi="", runopt="gen"):
        """
        Reference:
            https://larrucea.eu/compute-ir-raman-spectra-qe/

        General procedure of calculation IR and Raman using ph.x mainly
            1. Optimize the wavefunction by performing an Self Consistent Field (scf) calculation with pw.x
            2. Calculate the vibrational frequencies (normal modes/phonons) with ph.x
            3. Extract the phonon information from ph.x output using dynmat.x
            4. Parse the dynmat.x output section that contains the spectra data (frequencies and intensities) and plot it with gnuplot, producing these two spectra:
        """
        self.phx_qmesh(mpi=mpi, runopt=runopt)
        self.dynmat(mpi=mpi, runopt=runopt)

    def fermi_surface(self, directory="tmp-qe-static", inpname="fermi-surface.in", output="fermi-surface.out", mpi="", runopt="gen"):
        """
        scf->nscf(with denser k points)->fs.x
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("dynmat.x calculation:\n")
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
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s fs.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def pp(self, directory="tmp-qe-static", prefix="pp", plot_num=0, iflag=3, output_format=5, mpi="", runopt="gen"):
        """
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
            with open(os.path.join(directory, prefix+"-"+table[plot_num]+".in"), 'w') as fout:
                self.pp_inputpp(fout, plot_num=plot_num, filplot=table[plot_num]+".dat")
                self.pp_plot(fout, output_format=output_format, iflag=iflag, filepp=table[plot_num]+".dat")
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=prefix+"-"+table[plot_num]+".in", output=prefix+"-"+table[plot_num]+".out", cmd="pp.x")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s pp.x < %s | tee %s" % (mpi, prefix+"-"+table[plot_num]+".in", prefix+"-"+table[plot_num]+".out"))
            os.chdir("../")


    def pp_inputpp(self, fout, plot_num, filplot):
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

    def pp_plot(self, fout, filepp, iflag=3, output_format=5, 
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
            fout.write("nx = %d\n", nx)
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
            fout.write("nx = %d, ny = %d\n" (nx, ny))
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

    def xspectra(self, directory="tmp-qe-static", inpname="xspectra.in", output="xspectra.out", mpi="", runopt="gen"):
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

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s xspectra.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")
    #

    def gen_yh(self,inpname, output, directory="tmp-qe-static", cmd="pw.x"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % (cmd, inpname, output))

