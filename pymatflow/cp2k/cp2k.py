"""
Overall abstraction of CP2K
"""
import numpy as np
import sys
import os
import shutil
from pymatflow.remote.server import server_handle

from pymatflow.cp2k.base.atom import cp2k_atom
from pymatflow.cp2k.base.debug import cp2k_debug
from pymatflow.cp2k.base.ext_restart import cp2k_ext_restart
from pymatflow.cp2k.base.farming import cp2k_farming
from pymatflow.cp2k.base.force_eval import cp2k_force_eval
from pymatflow.cp2k.base.glob import cp2k_glob
from pymatflow.cp2k.base.motion import cp2k_motion
from pymatflow.cp2k.base.multiple_force_evals import cp2k_multiple_force_evals
from pymatflow.cp2k.base.negf import cp2k_negf
from pymatflow.cp2k.base.optimize_basis import cp2k_optimize_basis
from pymatflow.cp2k.base.optimize_input import cp2k_optimize_input
from pymatflow.cp2k.base.swarm import cp2k_swarm
from pymatflow.cp2k.base.test import cp2k_test
from pymatflow.cp2k.base.vibrational_analysis import cp2k_vibrational_analysis


"""
"""

class cp2k:
    """
    Philosophy:
        I have to say the implementation way of cp2k.base.xxx and cp2k is actually
        not efficient in running.
        and we know that for a specific type of calculation like statis scf
        we only need GLOBAL and FORCE_EVAL. but we build static_run class
        inheriting from class cp2k, which will make it holds other modules
        like motion, farming, atom, multipole_force_evals, etc. yes that may
        be tedious from the perspective of view of coding. and in the before
        there is no class cp2k, every type of calculation if customly build
        only including the needed input module. however, I think to build a
        class like cp2k will providing more convinence for managing codes.
        And most importantly the code will run not very computationally
        consuming. So I consist on this.
    """
    def __init__(self):
        """
        """
        self.atom = cp2k_atom()
        self.debug = cp2k_debug()
        self.ext_restart = cp2k_ext_restart()
        self.farming = cp2k_farming()
        self.force_eval = cp2k_force_eval()
        self.glob = cp2k_glob()
        self.motion = cp2k_motion()
        self.multipole_force_evals = cp2k_multiple_force_evals()
        self.negf = cp2k_negf()
        self.optimize_basis = cp2k_optimize_basis()
        self.optimize_input = cp2k_optimize_input()
        self.swarm = cp2k_swarm()
        self.test = cp2k_test()
        self.vibrational_analysis = cp2k_vibrational_analysis()

        self._initialize()

    def _initialize(self):
        """ initialize the current object, do some default setting
        """
        self.run_params = {}
        self.set_run()

    def get_xyz(self, xyzfile):
        """
        :param xyzfile:
            a modified xyz formatted file(the second line specifies the cell of the
            system).
        """
        self.force_eval.subsys.xyz.get_xyz(xyzfile)

    def set_params(self, params={}):
        """
        Note:
            we should always use this function to set params in cp2k

            every item in params begin with a "XXX-" where XXX is the first level
            input section of cp2k, like "FORCE_EVAL", "ATOM", "MOTION".

            currently only FORCE_EVAL, MOTION, ATOM, VIBRATIONAL_ANALYSIS are reliable.
        """
        for item in params:
            if item.split("-")[0].upper() == "ATOM":
                self.atom.set_params({item.replace("ATOM-", ""): params[item]})
            elif item.split("-")[0].upper() == "DEBUG":
                self.debug.set_params({item: params[item]})
            elif item.split("-")[0].upper() == "EXT_RESTART":
                self.ext_restart.set_params({item: params[item]})
            elif item.split("-")[0].upper() == "FARMING":
                self.farming.set_params({item: params[item]})
            elif item.split("-")[0].upper() == "FORCE_EVAL":
                self.force_eval.set_params({item.replace("FORCE_EVAL-", ""): params[item]})
            elif item.split("-")[0].upper() == "GLOBAL":
                self.glob.set_params({item: params[item]})
            elif item.split("-")[0].upper() == "MOTION":
                self.motion.set_params({item.replace("MOTION-", ""): params[item]})
            elif item.split("-")[0].upper() == "MULTIPOLE_FORCE_EVALS":
                self.multipole_force_evals.set_params({item: params[item]})
            elif item.split("-")[0].upper() == "NEGF":
                self.negf.set_params({item: params[item]})
            elif item.split("-")[0].upper() == "OPTIMIZE_BASIS":
                self.optimize_basis.set_params({item: params[item]})
            elif item.split("-")[0].upper() == "OPTIMIZE_INPUT":
                self.optimize_input.set_params({item: params[item]})
            elif item.split("-")[0].upper() == "SWARM":
                self.swarm.set_params({item: params[item]})
            elif item.split("-")[0].upper() == "TEST":
                self.test.set_params({item: params[item]})
            elif item.split("-")[0].upper() == "VIBRATIONAL_ANALYSIS":
                self.vibrational_analysis.set_params({item.replace("VIBRATIONAL_ANALYSIS-", ""): params[item]})

    def set_vdw(self, usevdw=False):
        if usevdw == True:
            self.force_eval.dft.xc.vdw_potential.status = True
        else:
            self.force_eval.dft.xc.vdw_potential.status = False

    def set_printout(self, option=[]):
        """
        Note:
            responsible for the parseing of the printout_option
        :param option:
            1: printout pdos
            2: printout band
            3: printout electron densities
            4: printout electron local function(ELF)
            5: printout molecular orbitals
            6: printout molecular orbital cube files
            7: printout mulliken populaltion analysis
            8: printout cubes for generation of STM images
            9: printout cube file with total density(electrons+atomic core)
           10: printout v_hartree_cube
           11: printout v_xc_cube
           12: printout xray_diffraction_spectrum
           13: request a RESP fit of charges.
        """
        self.force_eval.dft.printout.status = True
        self.force_eval.properties.status = True

        if 1 in option:
            self.force_eval.dft.printout.pdos.status = True
        if 2 in option:
            self.force_eval.dft.printout.band_structure.status = True
            # simply set status to True will not work !!!!
            # you have to also set kpath throught the following commented function
            #self.force_eval.dft.printout.band_structure.set_band(kpath=kpath)
        if 3 in option:
            self.force_eval.dft.printout.e_density_cube.status = True
        if 4 in option:
            self.force_eval.dft.printout.elf_cube.status = True
        if 5 in option:
            self.force_eval.dft.printout.mo.status = True
        if 6 in option:
            self.force_eval.dft.printout.mo_cubes.status = True
        if 7 in option:
            self.force_eval.dft.printout.mulliken.status = True
        if 8 in option:
            self.force_eval.dft.printout.stm.status = True
        if 9 in option:
            self.force_eval.dft.printout.tot_density_cube.status = True
        if 10 in option:
            self.force_eval.dft.printout.v_hartree_cube.status = True
        if 11 in option:
            self.force_eval.dft.printout.v_xc_cube.status = True
        if 12 in option:
            self.force_eval.dft.printout.xray_diffraction_spectrum.status = True
        if 13 in option:
            self.force_eval.properties.resp.status = True

    def set_run(self, mpi="", server="pbs", jobname="cp2k", nodes=1, ppn=32):
        """ used to set  the parameters controlling the running of the task
        :param mpi: you can specify the mpi command here, it only has effect on native running

        """
        self.run_params["server"] = server
        self.run_params["mpi"] = mpi
        self.run_params["jobname"] = jobname
        self.run_params["nodes"] = nodes
        self.run_params["ppn"] = ppn
        #self.run_params["inpname"] = inpname
        #self.run_params["output"] = output


    def gen_yh(self, inpname, output, directory, cmd="cp2k.psmp"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".inp")[0]+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s -in %s | tee %s\n" % (cmd, inpname, output))


    def gen_pbs(self, inpname, output, directory, cmd="cp2k.popt", jobname="cp2k", nodes=1, ppn=32):
        """
        generating pbs job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".inp")[0]+".pbs"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % jobname)
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s -in %s > %s\n" % (cmd, inpname, output))
