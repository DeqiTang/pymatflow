#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from emuhelper.qe.base.control import qe_control
from emuhelper.qe.base.system import qe_system
from emuhelper.qe.base.electrons import qe_electrons
from emuhelper.qe.base.ions import qe_ions
from emuhelper.qe.base.arts import qe_arts


class static_run:
    """
    """
    def __init__(self, xyz_f):
        self.control = qe_control()
        self.system = qe_system()
        self.electrons = qe_electrons()
        self.arts = qe_arts(xyz_f)

        self.control.basic_setting("scf") 
        
        self.system.basic_setting(self.arts)

    def gen_input(self, directory="tmp-qe-static", inpname="static-scf.in"):
        """
        directory: a place for all the generated files
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        
        os.system("cp *.UPF %s/" % directory)

        with open(os.path.join(directory, inpname), 'w') as fout:
            self.control.to_in(fout)
            self.system.to_in(fout)
            self.electrons.to_in(fout)
            self.arts.to_in(fout)
    
    def run_scf(self, directory="tmp-qe-static", inpname="static-scf.in", output="static-scf.out"):
        os.chdir(directory)
        os.system("pw.x < %s | tee %s" % (inpname, output))
        os.chdir("../")

    def run_nscf(self, directory="tmp-qe-static", inpname="static-nscf.in", output="static-nscf.out"):
        os.chdir(directory)
        os.system("pw.x < %s | tee %s" % (inpname, output))
        os.chdir("../")
    
    def converge_ecutwfc(self, emin, emax, step, directory="tmp-qe-ecutwfc"):
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        os.system("cp *.UPF %s/" % directory)
    
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
        # run the simulation
        for i in range(n_test + 1):
            ecut_wfc = int(emin + i * step)
            inp_name = "ecutwfc-%d.in" % ecut_wfc
            out_f_name = "ecutwfc-%d.out" % ecut_wfc
            os.system("pw.x < %s | tee %s" % (inp_name, out_f_name))

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

        plt.plot(ecut_wfc_all, energy_all)
        plt.show()
        os.chdir("../")

        
    def converge_ecutrho(self, emin, emax, step, ecutwfc, directory="tmp-qe-ecutrho"):
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
        # run the simulation
        for i in range(n_test + 1):
            ecut_rho = int(emin + i * step)
            inp_name = "ecutrho-%d.in" % ecut_rho
            out_f_name = "ecutrho-%d.out" % ecut_rho
            os.system("pw.x < %s | tee %s" % (inp_name, out_f_name))
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

        plt.plot(ecut_rho_all, energy_all)
        plt.show()
        os.chdir("../")
    #
    def converge_kpoints(self,nk_min, nk_max, step=1, directory="tmp-qe-kpoints", ecutwfc=50):
        """
        test the energy convergenc against k-points

        currently only support automatic schme of K_POINTS
        and only nk1 = nk2 = nk3 are supported
        """
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
            self.system.params['ecutwfc'] = ecutwfc # use the previously converged ecutwfc
            self.arts.set_kpoints([nk, nk, nk, 0, 0, 0])
            with open(inp_name, 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)
        # run the simulation
        for i in range(n_test + 1):
            nk = nk_min + i * step
            inp_name = "kpoints-%d.in" % nk
            out_f_name = "kpoints-%d.out" % nk
            os.system("pw.x < %s | tee %s" % (inp_name, out_f_name))

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

        plt.plot(nk_all, energy_all)
        plt.show()
        os.chdir("../")  

    def converge_degauss(self,degauss_min, degauss_max, step, smearing='gauss', directory="tmp-qe-degauss", ecutwfc=50):
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
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        os.system("cp *.UPF %s/" % directory)
    
        os.chdir(directory)
        
        n_test = int((degauss_max - degauss_min) / step)
        for i in range(n_test + 1):
            degauss = degauss_min + i * step
            inp_name = "degauss-%d.in" % degauss
            self.control.params['outdir'] = './tmp-' + str(degauss)
            self.system.params['ecutwfc'] = ecutwfc # use the previously converged ecutwfc
            #self.arts.set_kpoints([nk, nk, nk, 0, 0, 0]) # use the previously convered kpoints(automatic)
            self.system.params['smearing'] = smearing
            self.system.params['degauss'] = degauss
            with open(inp_name, 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)
        # run the simulation
        for i in range(n_test + 1):
            degauss = degauss_min + i * step
            inp_name = "degauss-%d.in" % degauss
            out_f_name = "degauss-%d.out" % degauss
            os.system("pw.x < %s | tee %s" % (inp_name, out_f_name))

        # analyse the result
        for i in range(n_test + 1):
            degauss = degauss_min + i * step
            out_f_name = "degauss-%d.out" % degauss
            os.system("cat %s | grep '!    total energy' >> energy-degauss.data" % out_f_name)

        degauss_all = [degauss_min + i * step for i in range(n_test + 1)]
        energy_all = []
        with open("energy-degauss.data", 'r') as fin:
            for line in fin:
                energy_all.append(float(line.split()[4]))

        plt.plot(degauss_all, energy_all)
        plt.show()
        os.chdir("../")  

    def scf(self, directory="tmp-qe-static"):
        self.control.calculation("scf")

    def nscf(self, directory="tmp-qe-static", inpname="static-nscf.in"):
        """
        first check whether there is a previous scf running
        """
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("non-scf calculation:\n")
            print("  directory of previous scf calculattion not found!\n")
            sys.exit(1)

        self.control.calculation("nscf")
        self.arts.set_kpoints([6, 6, 6, 0, 0, 0])
        with open(os.path.join(directory, inpname), 'w') as fout:
            self.control.to_in(fout)
            self.system.to_in(fout)
            self.electrons.to_in(fout)
            self.arts.to_in(fout)
    
    def dos(self, directory="tmp-qe-static", inpname="static-dos.in", output="static-dos.out"):
        """
        first check whether there is a previous scf running
        """
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("dos calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        with open(os.path.join(directory, inpname), 'w') as fout:
            fout.write("&DOS\n")
            #fout.write("prefix = '%s'\n" % "pwscf")
            fout.write("outdir = '%s'\n" % self.control.params["outdir"])
            fout.write("fildos = 'output.dos'\n")
            fout.write("/\n")
            fout.write("\n")
        os.chdir(directory)
        os.system("dos.x < %s | tee %s" % (inpname, output))
        os.chdir("../")

    def bands(self, directory="tmp-qe-static", inpname1="static-bands.in", output1="static-bands.out", inpname2="bands.in", output2="bands.out"):
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

        self.control.calculation('bands')
        self.arts.set_kpoints(option="crystal_b")
        #self.arts.set_kpoints([6, 6, 6, 0, 0, 0])
        with open(os.path.join(directory, inpname1), 'w') as fout:
            self.control.to_in(fout)
            self.system.to_in(fout)
            self.electrons.to_in(fout)
            self.arts.to_in(fout)

        os.chdir(directory)
        os.system("pw.x < %s | tee %s" % (inpname1, output1))
        os.chdir("../")

        with open(os.path.join(directory, inpname2), 'w') as fout:
            fout.write("&BANDS\n")
            fout.write("outdir = '%s'\n" % ("./tmp"))
            fout.write("filband = '%s'\n" % ("bands.dat"))
            fout.write("lsym = %s\n" % (".true."))
            fout.write("/\n")
            fout.write("\n")

        os.chdir(directory)
        os.system("bands.x < %s | tee %s" % (inpname2, output2))
        os.chdir("../")
        

    def projwfc(self, directory="tmp-qe-static", inpname="static-projwfc.in", output="static-projwfc.out"):
        """
        &projwfc can using projwfc.x to calculate Lowdin charges, spilling 
        parameter, projected DOS
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("projwfc calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        with open(os.path.join(directory, inpname), 'w') as fout:
            fout.write("&PROJWFC\n")
            fout.write("outdir = '%s'\n" % self.control.params["outdir"])
            fout.write("filpdos = 'output.projwfc'\n")
            fout.write("/\n")
            fout.write("\n")
        os.chdir(directory)
        os.system("projwfc.x < %s | tee %s" % (inpname, output))
        os.chdir("../")

    def epsilon(self, directory="tmp-qe-static", inpname="epsilon.in", output="epsilon.out"):
        """
        References:
            https://gitlab.com/QEF/material-for-ljubljana-qe-summer-school/blob/master/Day-3/handson-day3-TDDFPT.pdf

        epsilon.x:
            calculation of absorption spectra in IPA(Independent Particle 
            Approximation).
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("epsilon.x calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        with open(os.path.join(directory, inpname), 'w') as fout:
            fout.write("&INPUTPP\n")
            fout.write("calculation = 'eps'\n")
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
        os.chdir(directory)
        os.system("epsilon.x < %s | tee %s" % (inpname, output))
        os.chdir("../")

    def turbo_davidson(self, directory="tmp-qe-static", inpname1="turbo-davidson.in", output1="turbo-davidson.out",
            inpname2="turbo-spectrum-davidson.in", output2="turbo-spectrum-davidson.out"):
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
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("turbo_davidson calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        with open(os.path.join(directory, inpname1), 'w') as fout:
            fout.write("&lr_input\n")
            fout.write("outdir = '%s'\n" % self.control.params["outdir"])
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
        os.chdir(directory)
        os.system("turbo_davidson.x < %s | tee %s" % (inpname1, output1))
        os.chdir("../")
    
        with open(os.path.join(directory, inpname2), 'w') as fout:
            fout.write("&lr_input\n")
            fout.write("outdir = %s\n" % self.control.params["outdir"])
            fout.write("td = 'davidson'\n")
            fout.write("epsil = 0.004\n")
            fout.write("start = 0.0d0\n")
            fout.write("end = 1.0d0\n")
            fout.write("increment = 0.0001d0\n")
            fout.write("eign_file = 'pwscf.eigen'\n")
            fout.write("/\n")
            fout.write("\n")
        os.chdir(directory)
        os.system("turbo_spectrum.x < %s | tee %s" % (inpname2, output2))
        os.chdir("../")
    
    def turbo_lanczos(self, directory="tmp-qe-static", inpname1="turbo-lanczos.in", output1="turbo-lanczos.out",
            inpname2="turbo-spectrum-lanczos.in", output2="turbo-spectrum-lanczos.out"):
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
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("turbo_lanczos calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        with open(os.path.join(directory, inpname1), 'w') as fout:
            fout.write("&lr_input\n")
            fout.write("outdir = '%s'\n" % self.control.params["outdir"])
            fout.write("restart_step = 100\n")
            fout.write("restart = .false.\n")
            fout.write("/\n")
            fout.write("\n")
            fout.write("&lr_control\n")
            fout.write("intermax = 1500\n")
            fout.write("ipol = 1\n")
            fout.write("/\n")
            fout.write("\n")
        os.chdir(directory)
        os.system("turbo_lanczos.x < %s | tee %s" % (inpname1, output1))
        os.chdir("../")
    
        with open(os.path.join(directory, inpname2), 'w') as fout:
            fout.write("&lr_input\n")
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
        os.chdir(directory)
        os.system("turbo_spectrum.x < %s | tee %s" % (inpname2, output2))
        os.chdir("../")

    def phx_qmesh(self, directory="tmp-qe-static", inpname="phx-qmesh.in", output="phx-qmesh.out", dynamat_file="phx-qmesh.dyn"):
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
        with open(os.path.join(directory, inpname), 'w') as fout:
            fout.write("&inputph\n")
            fout.write("tr2_ph = 1.0d-14\n")
            fout.write("ldisp = .true.\n")
            fout.write("nq1 = 2 \n") # 4
            fout.write("nq2 = 2 \n") # 4
            fout.write("nq3 = 2 \n") # 4
            fout.write("outdir = './tmp'\n")
            fout.write("fildyn = '%s'\n" % dynamat_file)
            fout.write("/\n")
            fout.write("\n")
        os.chdir(directory)
        os.system("ph.x < %s | tee %s" % (inpname, output))
        os.chdir("../")

    def q2r(self, directory="tmp-qe-static", inpname="q2r.in", output="q2r.out", dynamat_file="phx-qmesh.dyn"):
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
        with open(os.path.join(directory, inpname), 'w') as fout:
            fout.write("&input\n")
            fout.write("fildyn = '%s'\n" % dynamat_file) # Dynamical matrices from the phonon calculation
            fout.write("zasr = 'simple'\n") # A way to impose the acoustic sum rule
            fout.write("flfrc = ifc.fc\n") # Output file of the interatomic force constants
            fout.write("/\n")
            fout.write("\n")
        os.chdir(directory)
        os.system("q2r.x < %s | tee %s" % (inpname, output))
        os.chdir("../")

    def matdyn(self, directory="tmp-qe-static", inpname="matdyn.in", output="matdyn.out", ifc_file="ifc.fc"):
        """
        matdyn.x
            Calculate phonons at generic q points using IFC
        """
        nqpoints = 2
        qpoints = [
                [0.0, 0.0, 0.0],
                [0.012658, 0.0, 0.012658]
                ]
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("matdyn.x calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        with open(os.path.join(directory, inpname), 'w') as fout:
            fout.write("&input\n")
            fout.write("asr = 'simple'\n") # Acoustic sum rule
            fout.write("flfrc = %s\n" % ifc_file) # File with IFC's
            fout.write("flfrq = 'frequencies.freq'\n") # Output file with the frequencies
            fout.write("/\n")
            fout.write("%d\n" % nqpoints) # Number of q points
            for i in range(nqpoints):
                fout.write("%f %f %f\n" % (qpoints[i][0], qpoints[i][1], qpoints[i][2]))
        os.chdir(directory)
        os.system("matdyn.x < %s | tee %s" % (inpname, output))
        os.chdir("../")

    def plotband(self, directory="tmp-qe-static", inpname="plotband.in", output="plotband.out"):
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
        with open(os.path.join(directory, inpname), 'w') as fout:
            fout.write("%s\n" % frequencies_file) # Input file with the frequencies at various q
            fout.write("%f %f\n" % (freq_range_min, freq_range_max)) # Range of frequencies for a visualization
            fout.write("freq.plot\n") # Output file with frequencies which will be used for plot
            fout.write("freq.ps\n") # Plot of the dispersion
            fout.write("0.0\n") # Fermi level (needed only for band structure plot)
            fout.write("100.0 0.0\n") # Freq. step and reference freq. on the plot freq.ps
        os.chdir(directory)
        os.system("plotband.x < %s | tee %s" % (inpname, output))
        os.chdir("../")


    def phx_gamma(self, directory="tmp-qe-static", inpname="phx-gamma.in", output="phx-gamma.out", dynamat_file="phx-gamma.dyn"):
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
        with open(os.path.join(directory, inpname), 'w') as fout:
            fout.write("&inputph\n")
            fout.write("tr2_ph = 1.0d-14\n")
            fout.write("oudir = './tmp'")
            fout.write("fildyn = '%s'\n" % dynamat_file)
            fout.write("/\n")
            fout.write("0.0 0.0 0.0\n")
            fout.write("\n")
        os.chdir(directory)
        os.system("ph.x < %s | tee %s" % (inpname, output))
        os.chdir("../")

    def dynmat(self, directory="tmp-qe-static", inpname="dynmat.in", output="dynmat.out", dynamat_file="phx-qmesh.dyn"):
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
        with open(os.path.join(directory, inpname), 'w') as fout:
            fout.write("&input\n")
            fout.write("fildyn = '%s'\n" % dynamat_file) # File containing the dynamical matrix
            fout.write("asr = 'simple'\n")
            fout.write("/\n")
            fout.write("\n")
        os.chdir(directory)
        os.system("dynmat.x < %s | tee %s" % (inpname, output))
        os.chdir("../")

    def ir_raman(self, directory="tmp-qe-static"):
        """
        Reference:
            https://larrucea.eu/compute-ir-raman-spectra-qe/

        General procedure of calculation IR and Raman using ph.x mainly
            1. Optimize the wavefunction by performing an Self Consistent Field (scf) calculation with pw.x
            2. Calculate the vibrational frequencies (normal modes/phonons) with ph.x
            3. Extract the phonon information from ph.x output using dynmat.x
            4. Parse the dynmat.x output section that contains the spectra data (frequencies and intensities) and plot it with gnuplot, producing these two spectra:
        """
        self.phx_qmesh()
        self.dynmat()
