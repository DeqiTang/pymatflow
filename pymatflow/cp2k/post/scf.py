# ==============================================================================
import datetime
import matplotlib.pyplot as plt

from pymatflow.base.atom import Atom
from pymatflow.base.xyz import base_xyz

class scf_out:
    """
    """
    def __init__(self):
        """
        output is the output file of scf run
        """
        self.file = None
        self.scf_params = {}
        self.run_info = {}

    def get_info(self, file):
        """
        get the general information of scf run from scf run output file
        which is now stored in self.lines
        """
        self.clean()

        self.file = file
        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()
        self.get_scf_params_and_run_info()

    def clean(self):
        self.file = None
        self.scf_params = {}
        self.run_info = {}

    #
    def get_scf_params_and_run_info(self):
        """
        self.run_info[]
            start_time: the task start time
            stop_time: the task stop time
            scf_energies: all the energies during the scf procedure
            #fermi_energy: fermi energy of the system (if output)

        """
        self.run_info["scf_energies"] = []

        for i in range(len(self.lines)):
            # if it is an empty line continue to next line
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0] == "****" and self.lines[i].split()[1] == "****" and self.lines[i].split()[5] == "STARTED":
                self.run_info["start_time"] = self.lines[i].split("\n")[0]
            elif self.lines[i].split()[0] == "****" and self.lines[i].split()[1] == "****" and self.lines[i].split()[5] == "ENDED":
                self.run_info["stop_time"] = self.lines[i].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Basis":
                self.scf_params["BASIS_SET_FILE_NAME"] = self.lines[i].split()[-1].split("\n")[0]
                self.scf_params["POTENTIAL_FILE_NAME"] = self.lines[i+1].split()[-1].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Project":
                self.scf_params["PROJECT_NAME"] = self.lines[i].split()[-1].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Run":
                self.scf_params["RUN_TYPE"] = self.lines[i].split()[-1].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Global":
                self.scf_params["PRINT_LEVEL"] = self.lines[i].split()[-1].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Total" and self.lines[i].split()[2] == "number":
                self.run_info["mpi_processes"] = int(self.lines[i].split()[-1].split("\n")[0])
                self.run_info["cpu_model"] = self.lines[i+3].split("name")[1].split("\n")[0]
            elif self.lines[i].split()[0] == "-" and self.lines[i].split()[1] == "Atoms:":
                self.run_info["n_atom"] = int(self.lines[i].split()[-1].split("\n")[0])
                self.run_info["n_shell"] = int(self.lines[i+2].split()[-1].split("\n")[0])
            elif self.lines[i].split()[0] == "max_scf:":
                self.scf_params["MAX_SCF"] = int(self.lines[i].split()[-1].split("\n")[0])
                self.scf_params["MAX_SCF_HISTORY"] = int(self.lines[i+1].split()[-1].split("\n")[0])
                self.scf_params["MAX_DIIS"] = int(self.lines[i+2].split()[-1].split("\n")[0])
            elif self.lines[i].split()[0] == "eps_scf:":
                self.scf_params["EPS_SCF"] = float(self.lines[i].split()[1])
                self.scf_params["EPS_SCF_HISTORY"] = float(self.lines[i+1].split()[1])
                self.scf_params["EPS_DIIS"] = float(self.lines[i+2].split()[1])
            elif self.lines[i].split()[0] == "Mixing" and self.lines[i].split()[1] == 'method:':
                self.scf_params["MIXING"] = self.lines[i].split()[2]
            elif self.lines[i].split()[0] == "added" and self.lines[i].split()[1] == "MOs":
                self.scf_params["ADDED_MOS"] = int(self.lines[i].split()[2])
            elif self.lines[i].split()[0] == "Number" and self.lines[i].split()[1] == "of" and self.lines[i].split()[2] == "electrons:":
                self.run_info["n_electrons"] = int(self.lines[i].split()[3])
                self.run_info["n_occcupied_orbital"] = int(self.lines[i+1].split()[4])
                self.run_info["n_molecular_orbital"] = int(self.lines[i+2].split()[4])
                self.run_info["n_orbital_function"] = int(self.lines[i+4].split()[4])
            elif self.lines[i].split()[0] == "Step" and self.lines[i].split()[1] == "Update" and self.lines[i].split()[-1] == "Change":
                self.run_info["scf_head_line"] = i
            elif self.lines[i].split()[0] == "***" and self.lines[i].split()[1] == "SCF" and self.lines[i].split()[2] == "run" and self.lines[i].split()[3] == "converged":
                self.run_info["scf_steps"] = int(self.lines[i].split()[5])
            elif self.lines[i].split()[0] ==  "Fermi" and self.lines[i].split()[1] == "energy:":
                self.run_info["fermi_energy"] = float(self.lines[i].split()[2])
            elif self.lines[i].split()[0] == "Total" and self.lines[i].split()[1] == "energy:":
                self.run_info["final_scf_energy"] = float(self.lines[i].split()[2])
            elif self.lines[i].split()[0] == "ENERGY|" and self.lines[i].split()[1] == "Total" and self.lines[i].split()[2] == "FORCE_EVAL":
                self.run_info["free_energy"] = float(self.lines[i].split()[8])
            else:
                pass
        # now we obtain the total energy for each dft scf step
        # note: cp2k will output the scf DFT energy for each scf step.
        # and the final_scf_energy is equal to the last value of scf_energies
        # the value of free_energy if list by ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):
        # which is the total DFT energy plus the Electronic entropic energy(very small).
        begin = self.run_info["scf_head_line"] + 2
        self.run_info["scf_energies"] = []
        for i in range(self.run_info["scf_steps"]):
            self.run_info["scf_energies"].append(float(self.lines[begin+i].split()[5]))
        # ----------------------------------------------------------------------
        # get the xyz structure from information extracted above:
        # WARNING: in low level print of cp2k, there is no structure coordinates
        # in the output file
        # ----------------------------------------------------------------------

    # def
    #
    
