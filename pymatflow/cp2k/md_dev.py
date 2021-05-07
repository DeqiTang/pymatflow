"""
Molecular Dynamics calculation
"""
import numpy as np
import sys
import os
import shutil


from pymatflow.remote.server import server_handle
from pymatflow.cp2k.cp2k_dev import Cp2k

"""
Reference:
    https://pc2.uni-paderborn.de/teaching/trainings/cp2k-tutorial/


TODO: implementing VIBRATIONAL SPECTRA calculating following this tutorial:
    https://pc2.uni-paderborn.de/teaching/trainings/cp2k-tutorial/
    it will use cp2k and Travis tools.
"""

class MdRun(Cp2k):
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
        super().__init__()

        self.set_params({
            "global-run_type": "MD",
            "force_eval-dft-mgrid-cutoff": 100,
            "force_eval-dft-mgrid-rel_cutoff": 60,
        })

        self.set_section_status({
            "motion-md": True,
        })


    def aimd(self, directory="tmp-cp2k-aimd", inpname="aimd.inp", output="aimd.out", runopt="gen", auto=0):
        """
        :param directory:
            directory is and path where the calculation will happen.
        :param inpname:
            input filename for the cp2k
        :param output:
            output filename for the cp2k
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.xyz.file, os.path.join(directory, os.path.basename(self.xyz.file)))

            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write(self.sections["global"].to_string())
                fout.write(self.sections["force_eval"].to_string())
                fout.write(self.sections["motion"].to_string())

            # gen server job comit file
            self.gen_llhpc(cmd="$PMF_CP2K", directory=directory, inpname=inpname, output=output)
            # gen pbs server job comit file
            self.gen_pbs(cmd="$PMF_CP2K", directory=directory, inpname=inpname, output=output, jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cdcloud server job comit file
            self.gen_cdcloud(cmd="$PMF_CP2K", directory=directory, inpname=inpname, output=output)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_CP2K -in %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="aimd", server=self.run_params["server"])

    def ir_spectra(self):
        """
        if you are calculating ir spectra, you have to
        have the dipole information for the molecule
        available in the simulated trajectory.
        that is realized by FORCE_EVAL%DFT%LOCALIZE

        Reference:
            http://www.travis-analyzer.de/
        """
        self.set_section_status({
            "force_eval-dft-localize": True,
        })


    def vib(self, directory="tmp-cp2k-md-vib", inpname="md.inp", output="md.out", mpi="", runopt="gen"):
        """
        :param directory:
            directory is and path where the calculation will happen.
        :param inpname:
            input filename for the cp2k
        :param output:
            output filename for the cp2k

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
            shutil.copyfile(self.xyz.file, os.path.join(directory, os.path.basename(self.xyz.file)))

            #self.force_eval.set_params(force_eval)
            #self.motion.set_params(motion)

            # --------------------------------------------
            # 1. trajectory -> a): massive equilibration    2000 steps
            # --------------------------------------------
            self.set_params({
                "FORCE_EVAL-DFT-QS-EPS_DEFAULT": 1.0e-12, # a good compromise
                "FORCE_EVAL-DFT-SCF-EPS_SCF": 1.0e-6,
                })
            # stop print wfn per md step to save time
            self.set_params({
                "FORCE_EVAL-DFT-SCF-PRINT-RESTART-EACH-MD": 0,
                })
            # Smoothing mitigates the break of translational invariance due to the plane waves
            # For cutoffs < 600 Ry (as we all use), this is absolutely mandatory
            self.set_section_status({
                "force_eval-dft-xc-xc_grid": True,
            })

            # setting of vde potential

            # set motion
            self.set_params({
                "MOTION-MD-ENSEMBLE": "NVT",
                "MOTION-MD-TIMESTEP": 0.5, # mandatory for some spectra calculation, like
                "MOTION-MD-STEPS": 200, #2000,
                "MOTION-MD-THERMOSTAT-TYPE": "NOSE",
                "MOTION-MD-THERMOSTAT-REGION": "MASSIVE",
                "MOTION-MD-THERMOSTAT-NOSE-TIMECON": 10.0,
                "MOTION-MD-TEMPERATURE": 300,
                })

            # Stop spamming all kinds of restart backup / history files
            # Only one single restart file, which is written in every MD step
            self.set_params({
                "MOTION-PRINT-RESTART-BACKUP_COPIES": 0,
                "MOTION-PRINT-RESTART-EACH-MD": 1,
                "MOTION-PRINT-RESTART_HISTORY-EACH-MD": 0,
                })
            with open(os.path.join(directory, "md-massive-equilibration.inp"), 'w') as fout:
                fout.write(self.sections["global"].to_string())
                fout.write(self.sections["force_eval"].to_string())
                fout.write(self.sections["motion"].to_string())
            # --------------------------------------------------------
            # 1. trajectory -> b): non-massive equilibration     20000 step
            # --------------------------------------------------------

            # restart from last WFN instead of initial guess
            self.set_params({
                "FORCE_EVAL-DFT-SCF-SCF_GUESS": "RESTART",
                })
            # Remove the MASSIVE after the first equilibration phase
            # Weaker thermostat coupling (time constant 100 fs) for second
            # equilibration and production run Strong thermostat coupling
            # might distort the dynamics and spectra...
            self.set_params({
                "MOTION-MD-ENSEMBLE": "NVT",
                "MOTION-MD-TIMESTEP": 0.5,
                "MOTION-MD-STEPS": 500,  # 20000,
                "MOTION-MD-THERMOSTAT-TYPE": "NOSE",
                "MOTION-MD-THERMOSTAT-REGION": None,
                "MOTION-MD-THERMOSTAT-NOSE-TIMECON": 100.0,
                })

            # need the EXT_RESTART block
            self..set_params({
                "EXT_RESTART-EXTERNAL_FILE": "ab-initio-1.restart",
                "EXT_RESTART-RESTART_THERMOSTAT": "FALSE",
                })
            with open(os.path.join(directory, "md-non-massive-equilibration.inp"), 'w') as fout:
                fout.write(self.sections["global"].to_string())
                fout.write(self.sections["force_eval"].to_string())
                fout.write(self.sections["motion"].to_string())
                fout.write(self.sections["ext_restart"].to_string())

            # ------------------------------------------------
            # 1. trajectory -> c): production run    60000 step
            # ------------------------------------------------

            self.set_params({
                "MOTION-MD-ENSEMBLE": "NVT",
                "MOTION-MD-TIMESTEP": 0.5,
                "MOTION-MD-STEPS": 1000, #60000,
                "MOTION-MD-THERMOSTAT-TYPE": "NOSE",
                "MOTION-MD-THERMOSTAT-REGION": None,
                "MOTION-MD-THERMOSTAT-NOSE-TIMECON": 100.0,
                })

            # need the EXT_RESTART block
            self.set_params({
                "EXT_RESTART-EXTERNAL_FILE": "ab-initio-1.restart",
                "EXT_RESTART-RESTART_THERMOSTAT": None,
                })
            with open(os.path.join(directory, "md-production-run.inp"), 'w') as fout:
                fout.write(self.sections["global"].to_string())
                fout.write(self.sections["force_eval"].to_string())
                fout.write(self.sections["motion"].to_string())
                fout.write(self.sections["ext_restart"].to_string())

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
            self.set_section_status({
                "force_eval-dft-periodic_efied": True,
            })

            self.set_params({
                "FORCE_EVAL-DFT-PERIODIC_EFIELD-INTENSITY": 5.0E-3,
                "FORCE_EVAL-DFT-PERIODIC_EFIELD-POLARIZATION": [1.0, 0.0, 0.0],
                })

            # Write the electron density in each MD step to a CUBE trajectory
            self.set_section_status({
                "force_eval-dft-print": True,
                "force_eval-dft-print-e_density_cube": True, # STRIDE 1 1 1 is vital for Voronoi integration
            }) 

            # Don‘t compute Wannier centers if you don‘t have to (can waste a lot of time if CRAZY does not converge)
            # Here we choose not to compute the Wannier centers
            self.set_section_status({
                "force_eval-dft-localize": False
            })

            # Follow the pre‐computed reference trajectory instead of doing a true MD
            # Make sure to specify the correct reference trajectory file name
            # Enter the FIRST_SNAPSHOT and STEPS according to your needs
            # EVAL_ENERGY _FORCES is important to re‐compute the electron structure
            self.set_params({
                "MOTION-MD-ENSEMBLE": "REFTRAJ",
                "MOTION-MD-STEPS": 1024,
                "MOTION-MD-TIMESTEP": None,
                "MOTION-MD-THERMOSTAT-TYPE": None,
                "MOTION-MD-THERMOSTAT-REGION": None,
                "MOTION-MD-THERMOSTAT-NOSE-TIMECON": None,
                "MOTION-MD-REFTRAJ-EVAL_ENERGY_FORCES": "TRUE",
                "MOTION-MD-REFTRAJ-FIRST_SNAPSHOT": 1,
                "MOTION-MD-REFTRAJ-TRAJ_FILE_NAME": "ab-initio-pos-1.xyz",
                })

            # This time: No restart files at all, because we just follow the reference trajectory
            self.set_params({
                "MOTION-PRINT-RESTART-BACKUP_COPIES": None,
                "MOTION-PRINT-RESTART-EACH-MD": 0,
                "MOTION-PRINT-RESTART_HISTORY-EACH-MD": 0,
                })

            with open(os.path.join(directory, "md-electron-density.inp"), 'w') as fout:
                fout.write(self.sections["global"].to_string())
                fout.write(self.sections["force_eval"].to_string())
                fout.write(self.sections["motion"].to_string())


            # -----------------------------------------------
            # 3.) Compress Volumetric Trajectories (Optional)
            # -----------------------------------------------


            # gen server job comit file

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("%s $PMF_CP2K -in %s | tee %s" % (mpi, inpname, output))
            os.system("%s $PMF_CP2K -in %s | tee %s" % (mpi, "md-massive-equilibration.inp", "md-massive-equilibration.out"))
            os.system("%s $PMF_CP2K -in %s | tee %s" % (mpi, "md-non-massive-equilibration.inp", "md-non-massive-equilibration.out"))
            os.system("%s $PMF_CP2K -in %s | tee %s" % (mpi, "md-production-run.inp", "md-production-run.out"))
            os.system("%s $PMF_CP2K -in %s | tee %s" % (mpi, "md-electron-density.inp", "md-electron-density.out"))
            os.chdir("../")

    #

    def metadynamics(self, directory="tmp-cp2k-metadynamics", inpname="metadynamics.inp", output="metadynamics.out", runopt="gen", auto=0):
        """
        :param directory:
            directory is and path where the calculation will happen.
        :param inpname:
            input filename for the cp2k
        :param output:
            output filename for the cp2k
        """
        self.set_section_status({
            "motion-free_energy": True,
            "motion-free_energy-metadyn": True,
            "motion-free_energy-metadyn-print": True,
        })
        
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write(self.sections["global"].to_string())
                fout.write(self.sections["force_eval"].to_string())
                fout.write(self.sections["motion"].to_string())

            # gen server job comit file
            self.gen_llhpc(cmd="$PMF_CP2K", directory=directory, inpname=inpname, output=output)
            # gen pbs server job comit file
            self.gen_pbs(cmd="$PMF_CP2K", directory=directory, inpname=inpname, output=output, jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cdcloud server job comit file
            self.gen_cdcloud(cmd="$PMF_CP2K", directory=directory, inpname=inpname, output=output)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_CP2K -in %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="metadynamics", server=self.run_params["server"])
