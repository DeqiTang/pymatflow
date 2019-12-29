#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.cp2k.base.glob import cp2k_glob
from pymatflow.cp2k.base.force_eval import cp2k_force_eval
from pymatflow.cp2k.base.motion import cp2k_motion
from pymatflow.cp2k.base.ext_restart import cp2k_ext_restart

"""
Reference:
    https://pc2.uni-paderborn.de/teaching/trainings/cp2k-tutorial/


TODO: implementing VIBRATIONAL SPECTRA calculating following this tutorial:
    https://pc2.uni-paderborn.de/teaching/trainings/cp2k-tutorial/
    it will use cp2k and Travis tools.
"""

class md_run:
    """
    Note:
        md_run is the class as an agent for Molecular Dynamics running. currently 
        implemented md type includes AIMD.
    TODO:
        implement QMMM and classic MD.
    """
    def __init__(self):
        """
        xyz_f:
            a modified xyz formatted file(the second line specifies the cell of the 
            system).
        """
        self.glob = cp2k_glob()
        self.force_eval = cp2k_force_eval()
        self.motion = cp2k_motion()
        self.ext_restart = cp2k_ext_restart()

        self.glob.basic_setting(run_type="MD")
        self.force_eval.basic_setting()

        self.motion.set_type("MD")

    def get_xyz(self, xyzfile):
        """
        xyz_f:
            a modified xyz formatted file(the second line specifies the cell of the 
            system).
        """
        self.force_eval.subsys.xyz.get_xyz(xyzfile)

    def set_params(self, force_eval={}, motion={}):
        """
        force_eval:
            allowing control of FORCE_EVAL/... parameters by user
        motion:
            allowing control of MOTION/... parameters by user
        """
        self.force_eval.set_params(force_eval)
        self.motion.set_params(motion)

    def aimd(self, directory="tmp-cp2k-aimd", inpname="aimd.inp", output="aimd.out", mpi="", runopt="gen"):
        """
        directory:
            directory is and path where the calculation will happen.
        inpname:
            input filename for the cp2k
        output:
            output filename for the cp2k
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.motion.to_input(fout)

            # gen server job comit file
            self.gen_yh(cmd="cp2k.popt", directory=directory, inpname=inpname, output=output)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s cp2k.psmp -in %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")   


    def ir_spectra(self):
        """
        if you are calculating ir spectra, you have to
        have the dipole information for the molecule 
        available in the simulated trajectory.
        that is realized by FORCE_EVAL%DFT%LOCALIZE

        Reference:
            http://www.travis-analyzer.de/
        """
        self.force_eval.dft.localize.status = True

    def vib(self, directory="tmp-cp2k-md-vib", inpname="md.inp", output="md.out", mpi="", runopt="gen"):
        """
        directory:
            directory is and path where the calculation will happen.
        inpname:
            input filename for the cp2k
        output:
            output filename for the cp2k
        force_eval:
            allowing control of FORCE_EVAL/... parameters by user
        motion:
            allowing control of MOTION/... parameters by user

        Note:
            standard: IR  | Raman
              chiral: VCD | ROA

            based on AIMD
            most IR/Raman spectra from AIMD use wannier centers. but it has 
            disadvantages: huge computational overhead, not guranteed to
            converge, only works for electrid dipole, cannot reproduce 
            quadrupole(required for ROA)
            the idea behind Travis is to completely drop wannier localization
            and integrate molecular dipole via Voronoi instead.

            usually it is hard to obtain VCD and ROA from AIMD, because
            both require the magnetic dipole moment.
            Travis has a solution.

        General workflow:
            cp2k    -> preparation
            cp2k    -> simulate trajectory
            cp2k    -> obtain electron density trajectories w/ ext. field
            bqbtool -> compress volumetric trajectories(optional)
            travis  -> solve current PDE, perform Voronoi integration
            travis  -> compute spectra from EMP property files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, self.force_eval.subsys.xyz.file))

            self.force_eval.set_params(force_eval)
            self.motion.set_params(motion)

            # --------------------------------------------
            # 1. trajectory -> a): massive equilibration    2000 steps
            # --------------------------------------------
            self.force_eval.set_params({
                "DFT-QS-EPS_DEFAULT": 1.0e-12, # a good compromise
                "DFT-SCF-EPS_SCF": 1.0e-6,
                })
            # stop print wfn per md step to save time
            self.force_eval.set_params({
                "DFT-SCF-PRINT-RESTART-EACH-MD": 0,
                })
            # Smoothing mitigates the break of translational invariance due to the plane waves
            # For cutoffs < 600 Ry (as we all use), this is absolutely mandatory
            self.force_eval.dft.xc.xc_grid.status = True

            # setting of vde potential

            # set motion
            self.motion.set_params({
                "MD-ENSEMBLE": "NVT",
                "MD-TIMESTEP": 0.5, # mandatory for some spectra calculation, like
                "MD-STEPS": 200, #2000,
                "MD-THERMOSTAT-TYPE": "NOSE",
                "MD-THERMOSTAT-REGION": "MASSIVE",
                "MD-THERMOSTAT-NOSE-TIMECON": 10.0,
                "MD-TEMPERATURE": 300,
                })
            
            # Stop spamming all kinds of restart backup / history files
            # Only one single restart file, which is written in every MD step
            self.motion.set_params({
                "PRINT-RESTART-BACKUP_COPIES": 0,
                "PRINT-RESTART-EACH-MD": 1,
                "PRINT-RESTART_HISTORY-EACH-MD": 0,
                })
            with open(os.path.join(directory, "md-massive-equilibration.inp"), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.motion.to_input(fout)
            # --------------------------------------------------------
            # 1. trajectory -> b): non-massive equilibration     20000 step
            # --------------------------------------------------------
            
            # restart from last WFN instead of initial guess
            self.force_eval.set_params({
                "DFT-SCF-SCF_GUESS": "RESTART",
                })
            # Remove the MASSIVE after the first equilibration phase
            # Weaker thermostat coupling (time constant 100 fs) for second 
            # equilibration and production run Strong thermostat coupling 
            # might distort the dynamics and spectra...
            self.motion.set_params({
                "MD-ENSEMBLE": "NVT",
                "MD-TIMESTEP": 0.5,
                "MD-STEPS": 500,  # 20000,
                "MD-THERMOSTAT-TYPE": "NOSE",
                "MD-THERMOSTAT-REGION": None,
                "MD-THERMOSTAT-NOSE-TIMECON": 100.0,
                })
            
            # need the EXT_RESTART block
            self.ext_restart.set_params({
                "EXTERNAL_FILE": "ab-initio-1.restart",
                "RESTART_THERMOSTAT": "FALSE",
                })
            with open(os.path.join(directory, "md-non-massive-equilibration.inp"), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.motion.to_input(fout)
                self.ext_restart.to_input(fout)

            # ------------------------------------------------
            # 1. trajectory -> c): production run    60000 step
            # ------------------------------------------------
            
            self.motion.set_params({
                "MD-ENSEMBLE": "NVT",
                "MD-TIMESTEP": 0.5,
                "MD-STEPS": 1000, #60000,
                "MD-THERMOSTAT-TYPE": "NOSE",
                "MD-THERMOSTAT-REGION": None,
                "MD-THERMOSTAT-NOSE-TIMECON": 100.0,
                })
            
            # need the EXT_RESTART block
            self.ext_restart.set_params({
                "EXTERNAL_FILE": "ab-initio-1.restart",
                "RESTART_THERMOSTAT": None,
                })
            with open(os.path.join(directory, "md-production-run.inp"), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.motion.to_input(fout)
                self.ext_restart.to_input(fout)
            # Now we have a 30 ps production trajectory which contains all motions for the spectra.


            # -------------------------------------------
            # 2.) Electron Density
            # -------------------------------------------

            # We traverse the trajectory again, and store volumetric electron density in every n‐th frame
            # For IR and VCD: Only one field‐free calculation required
            # For Raman and ROA: Field‐free + 3 field directions  4 runs
            # For IR and Raman: Sufficient to consider every 8 th frame(every 4.0 fs)
            # For VCD and ROA: Need every frame (every 0.5 fs)

            # need an external electric field which works with periodic systems
            # An external field strength of 5.0E‐3 a.u. is a good compromise 
            # between noise and linearity (corresponds to 2.5 * 10 9 V/m !)
            # POLARIZATION gives the field vector (here: positive X direction)
            self.force_eval.dft.periodic_efield.status = True

            self.force_eval.set_params({
                "DFT-PERIODIC_EFIELD-INTENSITY": 5.0E-3,
                "DFT-PERIODIC_EFIELD-POLARIZATION": [1.0, 0.0, 0.0],
                }) 

            # Write the electron density in each MD step to a CUBE trajectory
            self.force_eval.dft.printout.status = True
            self.force_eval.dft.printout.e_density_cube.status = True  # STRIDE 1 1 1 is vital for Voronoi integration

            # Don‘t compute Wannier centers if you don‘t have to (can waste a lot of time if CRAZY does not converge)
            # Here we choose not to compute the Wannier centers
            self.force_eval.dft.localize.status = False

            # Follow the pre‐computed reference trajectory instead of doing a true MD
            # Make sure to specify the correct reference trajectory file name
            # Enter the FIRST_SNAPSHOT and STEPS according to your needs
            # EVAL_ENERGY _FORCES is important to re‐compute the electron structure
            self.motion.set_params({
                "MD-ENSEMBLE": "REFTRAJ",
                "MD-STEPS": 1024,
                "MD-TIMESTEP": None,
                "MD-THERMOSTAT-TYPE": None,
                "MD-THERMOSTAT-REGION": None,
                "MD-THERMOSTAT-NOSE-TIMECON": None,
                "MD-REFTRAJ-EVAL_ENERGY_FORCES": "TRUE",
                "MD-REFTRAJ-FIRST_SNAPSHOT": 1,
                "MD-REFTRAJ-TRAJ_FILE_NAME": "ab-initio-pos-1.xyz",
                })

            # This time: No restart files at all, because we just follow the reference trajectory
            self.motion.set_params({
                "PRINT-RESTART-BACKUP_COPIES": None,
                "PRINT-RESTART-EACH-MD": 0,
                "PRINT-RESTART_HISTORY-EACH-MD": 0,
                })           

            with open(os.path.join(directory, "md-electron-density.inp"), 'w') as fout:
                self.glob.to_input(fout)
                self.force_eval.to_input(fout)
                self.motion.to_input(fout)


            # -----------------------------------------------
            # 3.) Compress Volumetric Trajectories (Optional)
            # -----------------------------------------------


            # gen server job comit file

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("%s cp2k.psmp -in %s | tee %s" % (mpi, inpname, output))
            os.system("%s cp2k.psmp -in %s | tee %s" % (mpi, "md-massive-equilibration.inp", "md-massive-equilibration.out"))
            os.system("%s cp2k.psmp -in %s | tee %s" % (mpi, "md-non-massive-equilibration.inp", "md-non-massive-equilibration.out"))
            os.system("%s cp2k.psmp -in %s | tee %s" % (mpi, "md-production-run.inp", "md-production-run.out"))
            os.system("%s cp2k.psmp -in %s | tee %s" % (mpi, "md-electron-density.inp", "md-electron-density.out"))
            os.chdir("../")
    


    def gen_yh(self,inpname, output, directory="tmp-cp2k-md", cmd="cp2k.psmp"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s -in %s | tee %s\n" % (cmd, inpname, output))
