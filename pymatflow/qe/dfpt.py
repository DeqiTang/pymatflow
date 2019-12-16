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


class dfpt_run:
    """
    About:
        dfpt_run implements the control over ph.x, dynmat.x, q2r.x,
        matdyn.x
        calculations based on them.
    Status:
        currently implemented calculation including:
            phx, q2r, matdyn, plotband,
            dynmat, ir_raman, fermi_surface,
    Note:
        ph.x calculation cannot start from pw.x data using Gamma-point
        tricks. so the static scf must be done not using Gamma kpoint
        scheme.
    
    occupations setting:
        
        ph.x calculation can not be going when the system is metallic, 

        but I found even when my fermi energy is in the gap, ph.x can sometimes fail to run. 
        this might result from use of smearing in the scf calculation. we should tray other type of occupation.

        sometimes ph.x can run when I use smearing type occupation,
        but somtimes it might warning it is metallic, and stop the calculation,
        even thought I found the fermi energy is actually in the gap(insulator).

        DFPT with the Blochl correction of occupation(tetrahedra) is not implemented
        but tetrahedra_opt and fixed is ok.

        for some systems if you do not use smearing occupations, during
        the scf ground state calculation qe will stop, signaling:
        'charge is wrong: smearing is needed'
        but actually in reality we know the system is an insulator.

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


    def phx(self, directory="tmp-qe-static", inpname="phx.in", output="phx.out", mpi="", runopt="gen",
            inputph={}, qpoints_option="gamma", qpoints=[2, 2, 2], polar="false"):
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

            PHonon: linear-response calculations(phonons, dielectric properties)
                (1) phonon frequencies and eigenvectors at a generic wave vector
                (2) dielectric tensor, effective charges, IR cross sections
                (3) interatomic force constants in real space
                (4) electron-phonon interaction coefficients for metals
                (5) nonresonant Raman cross sections
                (6) third-order anharmonic phonon lifetimes cross sections
        qpoints_option:
            'gamma' or 'qmesh'
        inputph:
            allowing setting of some self.inpuph parameters outside.
            be aware that parameters related to setting of qpoints
            and polar materials are not allowed to be set through
            inputph. namely nq1, nq2, nq3, and ldisp and epsil should
            not be in inputph
        """
        self.inputph = {
                "outdir": self.control.params["outdir"],
                "prefix": self.control.params["prefix"],
                "fildyn": "phx.dyn",
                "tr2_ph": 1.0e-14,
                "lrpa": None,
                "lnoloc": None,
                }
        # set the self.inputph through inputpp
        for forbid in ["nq1", "nq2", "nq3", "ldisp", "epsil"]:
            if forbid in inputph:
                print("=============================================\n")
                print("                 Warning!!!\n")
                print("=============================================\n")
                print("qe.dfpt.phx\n")
                print("the setting of [nq1, nq2, nq3, ldisp, epsil]\n")
                print("are not allowed through inputph\n")
        for item in inputph:
            self.inputph[item] = inputph[item]
        # ---------------------------------------------------------------
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("ph.x calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("ph.x calculation\n")
                fout.write("&inputph\n")
                for item in self.inputph:
                    if self.inputph[item] is not None:
                        if type(self.inputph[item]) == str:
                            fout.write("%s = '%s'\n" % (item, str(self.inputph[item])))
                        else:
                            fout.write("%s = %s\n" % (item, str(self.inputph[item])))
                # if polar = true will  set epsil = .true. 
                #to calculate and store the dielectric tensor and effective charge
                if polar.lower() == "yes":
                    fout.write("epsil = .true.\n")
                #
                if qpoints_option == "qmesh":
                    fout.write("nq1 = %d\n" % qpoints[0])
                    fout.write("nq2 = %d\n" % qpoints[1])
                    fout.write("nq3 = %d\n" % qpoints[2])
                    fout.write("ldisp = .true.\n") # must be set for q mesh calculation
                # get info about atom fixation
                nat_todo = 0
                for atom in self.arts.xyz.atoms:
                    if False in atom.fix:
                        nat_todo += 1
                        
                if nat_todo == self.arts.xyz.natom:
                    fout.write("nat_todo = 0\n") # displace all atoms
                    fout.write("nogg = .false.\n")
                else:
                    fout.write("nat_todo = %d\n" % nat_todo)
                    fout.write("nogg = .true.\n") 
                    # gamma_gamma tricks with nat_todo != 0 not available,
                    # so we must use nogg = .true.
                fout.write("/\n")
                if qpoints_option == "gamma":
                    fout.write("0.0 0.0 0.0\n")
                elif qpoints_option == "qmesh":
                    pass
                else:
                    print("=======================================\n")
                    print("              Warning !!!\n")
                    print("=======================================\n")
                    print("qe.dfpt.phx\n")
                    print("only qpoints_option = 'gamma' and 'qmesh'\n")
                    print("are supported\n")
                    sys.exit(1)
                # indicies of atom to be used in the calculation
                if nat_todo < self.arts.xyz.natom:
                    for i in range(self.arts.xyz.natom):
                        if False in self.arts.xyz.atoms[i].fix:
                            fout.write("%d " % (i+1))
                    fout.write("\n")
                # end incicies for atoms to be used in the calculation
                fout.write("\n")

            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="ph.x")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s ph.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def q2r(self, directory="tmp-qe-static", inpname="q2r.in", output="q2r.out", 
            dynamat_file="phx.dyn", ifc_file="q2r.fc", mpi="", runopt="gen", zasr='simple'):
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
            ifc_file="q2r.fc", mpi="", runopt="gen", asr='simple', 
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
                fout.write("flfrq = 'matdyn.freq'\n") # Output file with the frequencies
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

    def plotband(self, directory="tmp-qe-static", inpname="plotband.in", output="plotband.out", frequencies_file="matdyn.freq", mpi="", runopt="gen", freq_min=0, freq_max=600):
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
                fout.write("plotband-freq.plot\n") # Output file with frequencies which will be used for plot
                fout.write("plotband-freq.ps\n") # Plot of the dispersion
                fout.write("0.0\n") # Fermi level (needed only for band structure plot)
                fout.write("100.0 0.0\n") # Freq. step and reference freq. on the plot freq.ps
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="plotband.x")
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s plotband.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")


    def dynmat(self, directory="tmp-qe-static", inpname="dynmat.in", output="dynmat.out", mpi="", runopt="gen",
            fildyn="phx.dyn", asr="simple", qi=[0, 0, 0]):
        """
        imposing acoustic sum rule (ASR)
        extract the phonon information from ph.x output using dynmat.x(
        which can also be used to get IR and Raman.
        the generated fildyn.axsf fildyn.mold can be visualized by xcrysden
        and molden separately, and molden can visualize the vibration through
        fildyn.mold
        )
        Note:
            only used when the ph.x calculation was conducted using Gamma point but not the q mesh.
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
                fout.write("fildyn = '%s'\n" % fildyn) # File containing the dynamical matrix
                fout.write("asr = '%s'\n" % asr)
                fout.write("q(1) = %f\n" % qi[0])
                fout.write("q(2) = %f\n" % qi[1])
                fout.write("q(3) = %f\n" % qi[2])
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
        self.phx(mpi=mpi, runopt=runopt)
        self.dynmat(mpi=mpi, runopt=runopt)

    #
    def gen_yh(self, inpname, output, directory="tmp-qe-static", cmd="pw.x"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % (cmd, inpname, output))
