"""
TDDFPT calc
"""
import os
import sys
import shutil
import matplotlib.pyplot as plt


class tddfpt_run():
    """
    """
    def __init__(self):
        pass


    def set_epsilon(self, inputpp={}, energy_grid={}):
        self.inputpp_epsilon = {
                "calculation": 'eps',
                "prefix": "pwscf",
                "outdir": "./tmp",
                }
        self.energy_grid_epsilon = {
                "smeartype": 'gaussian',# type of smearing of the spectrum
                "intersmear": 0.1,# the valus of smearing in eV
                "wmin": 0.0,# minimum value of frequencies for a plot in eV
                "wmax": 15.0,# maximum value of frequencies for a plot in eV
                "nw": 1000,# number of points between wmin and wmax
                }
        for item in inputpp:
            self.inputpp_epsilon[item] = inputpp[item]
        for item in energy_grid:
            self.energy_grid_epsilon[item] = energy_grid[item]

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
                fout.write("calculation = '%s'\n" % self.inputpp_epsilon["calculation"])
                fout.write("prefix = '%s'\n" % self.inputpp_epsilon["prefix"])
                fout.write("outdir = '%s'\n" % self.inputpp_epsilon["outdir"])
                fout.write("/\n")
                fout.write("\n")
                fout.write("&ENERGY_GRID\n")
                fout.write("smeartype = '%s'\n" % self.energy_grid_epsilon["smeartype"]) # type of smearing of the spectrum
                fout.write("intersmear = %f\n" % self.energy_grid_epsilon["intersmear"])       # the valus of smearing in eV
                fout.write("wmin = %f\n" % self.energy_grid_epsilon["wmin"])             # minimum value of frequencies for a plot in eV
                fout.write("wmax = %f\n" % self.energy_grid_epsilon["wmax"])            # maximum value of frequencies for a plot in eV
                fout.write("nw = %f\n" % self.energy_grid_epsilon["nw"])              # number of points between wmin and wmax
                fout.write("/\n")
                fout.write("\n")

            # gen yh job submit script
            with open(os.path.join(directory, "epsilon.sub"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("yhrun -N 1 -n 24 epsilon.x < %s > %s\n" % (inpname, output))
            # gen pbs job submit script
            with open(os.path.join(directory, "epsilon.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % jobname)
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE epsilon.x < %s > %s\n" % (inpname, output))

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s epsilon.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")


    def set_turbo_davidson(self, lr_input={}, lr_dav={}):
        self.lr_input_td = {
                "outdir": "./tmp", #self.control.params["outdir"],
                "prefix": "pwscf", #self.control.params["prefix"],
                }
        self.lr_dav_td = {
                "if_dft_spectrum": False,
                "num_eign": 15,
                "num_init": 30,
                "num_basis_max": 90,
                "residue_conv_thr": 1.0e-6,
                "start": 0.0,
                "finish": 1.0,
                "step": 0.01,
                "broadening": 0.04,
                "reference": 0.3,
                }

        # set the self.lr_input_td through lr_input
        for item in lr_input:
            self.lr_input_td[item] = lr_input[item]

        # set the self.lr_dav_td through lr_dav
        for item in lr_dav:
            self.lr_dav_td[item] = lr_dav[item]

    def turbo_davidson(self, directory="tmp-qe-static", inpname1="turbo-davidson.in", output1="turbo-davidson.out",
            inpname2="turbo-spectrum-davidson.in", output2="turbo-spectrum-davidson.out", mpi="", runopt="gen",
            jobname="turbo-davidson", nodes=1, ppn=32):
        """
        references:
            https://gitlab.com/qef/material-for-ljubljana-qe-summer-school/blob/master/day-3/handson-day3-tddfpt.pdf

        turbo_davidson.x:
            calculate adsorption spectra of molecules using time-dependent
            density functinal perturbation theory(tddfpt).

            if if_dft_spectrum is set to .true. the result will be the same
            as calculated using epsilon.x, where there is no interaction.

            so set if_dft_spectrum to .false. to turn on the interaction.
            and you will get a shift of the peak compared to results calculated
            using epsilon.x or turbo_davidson.x(with if_dft_spectrum = .true.).

            when if_dft_spectrum is set to .true. turbo_davidson.x will produce
            a prefix-dft.eigen file, while a prefix.eigen file is produced if
            if_dft_spectrum is set to .false.

            we can also calculate absorption spectra using psudo-potential
            designed for b3lyp functional with turbo_davidson.x
            this way, we should set input_dft = 'b3lyp' in scf calcualtion
            and set d0psi_rs = .true. in input file for turbo_davidson.x

        turbo_spectrum.x:
            post-processing calculation of the spectrum
        note:
            turbotddft is not extended to metals, so we can only
            deal with insulators or semiconductors with turbo now.
            
            ltetra are not implemented now
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 warning !!!\n")
            print("===================================================\n")
            print("turbo_davidson calculation:\n")
            print("  directory of previous scf or nscf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname1), 'w') as fout:
                fout.write("&lr_input\n")
                fout.write("prefix = '%s'\n" % self.lr_input_td["prefix"])
                fout.write("outdir = '%s'\n" % self.lr_input_td["outdir"])
                #fout.write("wfcdir = ''")
                fout.write("/\n")
                fout.write("\n")
                fout.write("&lr_dav\n")
                fout.write("if_dft_spectrum = %s\n" % self.lr_dav_td["if_dft_spectrum"])
                fout.write("num_eign = %d\n" % self.lr_dav_td["num_eign"])
                fout.write("num_init = %d\n" % self.lr_dav_td["num_init"])
                fout.write("num_basis_max = %d\n" % self.lr_dav_td["num_basis_max"])
                fout.write("residue_conv_thr = %f\n" % self.lr_dav_td["residue_conv_thr"])
                fout.write("start = %f\n" % self.lr_dav_td["start"])
                fout.write("finish = %f\n" % self.lr_dav_td["finish"])
                fout.write("step = %f\n" % self.lr_dav_td["step"])
                fout.write("broadening = %f\n" % self.lr_dav_td["broadening"])
                fout.write("reference = %f\n" % self.lr_dav_td["reference"])
                fout.write("/\n")
                fout.write("\n")

            with open(os.path.join(directory, inpname2), 'w') as fout:
                fout.write("&lr_input\n")
                fout.write("prefix = '%s'\n" % self.lr_input_ts["prefix"])
                fout.write("outdir = %s\n" % self.lr_input_ts["outdir"])
                fout.write("td = 'davidson'\n")
                fout.write("epsil = %f\n" % self.lr_input_ts["epsil"])
                fout.write("start = %f\n" % self.lr_input_ts["start"])
                fout.write("end = %f\n" % self.lr_input_ts["end"])
                fout.write("increment = %f\n" % self.lr_input_ts["increment"])
                fout.write("eign_file = '%s'\n" % "pwscf.eigen")
                fout.write("/\n")
                fout.write("\n")
        
            # gen yh job submit script
            with open(os.path.join(directory, "turbo-davidson.sub"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("yhrun -N 1 -n 24 turbo_davidson.x < %s > %s\n" % (inpname1, output1))
                fout.write("yhrun -N 1 -n 24 turbo_spectrum.x < %s > %s\n" % (inpname2, output2))
            # gen pbs job submit script
            with open(os.path.join(directory, "turbo-davidson.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % jobname)
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE turbo_davidson.x < %s > %s\n" % (inpname1, output1))
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE turbo_spectrum.x < %s > %s\n" % (inpname2, output2))


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s turbo_davidson.x < %s | tee %s" % (mpi, inpname1, output1))
            os.system("%s turbo_spectrum.x < %s | tee %s" % (mpi, inpname2, output2))
            os.chdir("../")


    
    def set_turbo_spectrum(self, lr_input={}):
        self.lr_input_ts = {
                "outdir": "./tmp", #self.control.params["outdir"],
                "prefix": "pwscf", #self.control.params["prefix"],
                "restart_step": 100,
                "restart": False,
                #"td": "lanczos", # td not set by this function
                "epsil": 0.004, # The value of Lorenzian smearing in Ry
                "start": 0.0e0, # Minimum value of frequencies for a plot in Ry
                "end": 1.0e0, # Maximum value of frequencies for a plot in Ry
                "increment": 0.0001e0, # Frequency step in Ry
                "ipol": 1,  # Polarization direction (same as in turbo_lanczos.x)
                "eign_file": None, #self.control.params["prefix"]+".eigen",
                "intermax0": 1500, # number of calculated Lanczos coefficient
                "intermax": 20000, # number of extrapolated Lanczos coefficinet
                "extrapolation": 'osc', # Type of extrapolation (bi-constant)
                }

    def set_turbo_lanczos(self, lr_input={}, lr_control={}):
        self.lr_input_tl = {
                "outdir": "./tmp", #self.control.params["outdir"],
                "prefix": "pwscf", #self.control.params["prefix"],
                "restart_step": 100,
                "restart": False,
                }
        self.lr_control_tl = {
                "itermax": 1500,
                "ipol": 1,
                }

        # set the self.lr_input_tl through lr_input
        for item in lr_input:
            self.lr_input_tl[item] = lr_input[item]

        # set the self.lr_control_tl through lr_control
        for item in lr_control:
            self.lr_control_tl[item] = lr_control[item]

    def turbo_lanczos(self, directory="tmp-qe-static", inpname1="turbo-lanczos.in", output1="turbo-lanczos.out",
            inpname2="turbo-spectrum-lanczos.in", output2="turbo-spectrum-lanczos.out", mpi="", runopt="gen",
            jobname="turbo_lanczos", nodes=1, ppn=32):
        """
        references:
            https://gitlab.com/qef/material-for-ljubljana-qe-summer-school/blob/master/day-3/handson-day3-tddfpt.pdf

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
                fout.write("prefix = '%s'\n" % self.lr_input_tl["prefix"])
                fout.write("outdir = '%s'\n" % self.lr_input_tl["outdir"])
                fout.write("restart_step = %d\n" % self.lr_input_tl["restart_step"])
                fout.write("restart = %s\n" % self.lr_input_tl["restart"])
                fout.write("/\n")
                fout.write("\n")
                fout.write("&lr_control\n")
                fout.write("itermax = %d\n" % self.lr_control_tl["itermax"])
                fout.write("ipol = %d\n" % self.lr_control_tl["ipol"])
                fout.write("/\n")
                fout.write("\n")
        
            with open(os.path.join(directory, inpname2), 'w') as fout:
                fout.write("&lr_input\n")
                fout.write("prefix = '%s'\n" % self.lr_input_ts["prefix"])
                fout.write("outdir = %s\n" % self.lr_input_ts["outdir"])
                fout.write("intermax0 = %d\n" % self.lr_input_ts["intermax0"]) # Number of calculated Lanczos coefficient
                fout.write("intermax = %d\n" % self.lr_input_ts["intermax"]) # umber of extrapolated Lanczos coefficient
                fout.write("extrapolation = '%s'\n" % self.lr_input_ts["extrapolation"]) # Type of extrapolation (bi-constant)
                fout.write("epsil = %f\n" % self.lr_input_ts["epsil"]) # The value of Lorenzian smearing in Ry
                fout.write("start = %f\n" % self.lr_input_ts["start"]) # Minimum value of frequencies for a plot in Ry
                fout.write("end = %f\n" % self.lr_input_ts["end"]) # Maximum value of frequencies for a plot in Ry
                fout.write("increment = %f\n" % self.lr_input_ts["increment"]) # Frequency step in Ry
                fout.write("ipol = %d\n" % self.lr_input_ts["ipol"]) # Polarization direction (same as in turbo_lanczos.x)
                fout.write("/\n")
                fout.write("\n")

            # gen yh job submit script
            with open(os.path.join(directory, "turbo-lanczos.sub"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("yhrun -N 1 -n 24 turbo_lanczos.x < %s > %s\n" % (inpname1, output1))
                fout.write("yhrun -N 1 -n 24 turbo_spectrum.x < %s > %s\n" % (inpname2, output2))
            # gen pbs job submit script
            with open(os.path.join(directory, "turbo-lanczos.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % jobname)
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE turbo_lanczos.x < %s > %s\n" % (inpname1, output1))
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE turbo_spectrum.x < %s > %s\n" % (inpname2, output2))

        
        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s turbo_lanczos.x < %s | tee %s" % (mpi, inpname1, output1))
            os.system("%s turbo_spectrum.x < %s | tee %s" % (mpi, inpname2, output2))
            os.chdir("../")

