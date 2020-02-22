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
            if self.lines[i].split()[0] == "Program" and self.lines[i].split()[1] == "PWSCF" and self.lines[i].split()[3] == "starts":
                self.run_info["start_time"] = self.lines[i].split("\n")[0]
            elif self.lines[i].split()[0] == "This" and self.lines[i].split()[1] == "run" and self.lines[i].split()[3] == "terminated":
                self.run_info["stop_time"] = self.lines[i].split("\n")[0]
            elif self.lines[i].split()[0] == "Parallel" and self.lines[i].split()[-1] == "processors":
                self.run_info["processors"] = int(self.lines[i].split()[-2])
            elif self.lines[i].split()[0] == "MPI" and self.lines[i].split()[-1] == "nodes":
                self.run_info["nodes"] = int(self.lines[i].split()[-2])
            elif self.lines[i].split()[0] == "bravais-lattice" and self.lines[i].split()[1] == "index":
                self.scf_params["alat_au"] = float(self.lines[i+1].split()[4])
                self.scf_params["nat"] = int(self.lines[i+3].split()[4])
                self.scf_params["nelectron"] = float(self.lines[i+5].split()[4])
                self.scf_params["n_ks_state"] = int(self.lines[i+6].split("=")[1])
                self.scf_params["ecutwfc"] = int(float(self.lines[i+7].split()[3]))
                self.scf_params["ecutrho"] = int(float(self.lines[i+8].split()[4]))
                self.scf_params["conv_thr"] = float(self.lines[i+9].split()[3])
                self.scf_params["mixing_beta"] = float(self.lines[i+10].split()[3])
            elif self.lines[i].split()[0] == "crystal" and self.lines[i].split()[1] == "axes:" and self.lines[i].split()[-1] =="alat)":
                self.scf_params["cell_a_alat"] = []
                self.scf_params["cell_a_alat"].append([float(self.lines[i+1].split()[3]), float(self.lines[i+1].split()[4]), float(self.lines[i+1].split()[5])])
                self.scf_params["cell_a_alat"].append([float(self.lines[i+2].split()[3]), float(self.lines[i+2].split()[4]), float(self.lines[i+2].split()[5])])
                self.scf_params["cell_a_alat"].append([float(self.lines[i+3].split()[3]), float(self.lines[i+3].split()[4]), float(self.lines[i+3].split()[5])])
            elif self.lines[i].split()[0] == "reciprocal" and self.lines[i].split()[1] == "axes:" and self.lines[i].split()[-1] == "pi/alat)": # actually '2 pi/alat'
                self.scf_params["cell_b_2pi_alat"] = []
                self.scf_params["cell_b_2pi_alat"].append([float(self.lines[i+1].split()[3]), float(self.lines[i+1].split()[4]), float(self.lines[i+1].split()[5])])
                self.scf_params["cell_b_2pi_alat"].append([float(self.lines[i+2].split()[3]), float(self.lines[i+2].split()[4]), float(self.lines[i+2].split()[5])])
                self.scf_params["cell_b_2pi_alat"].append([float(self.lines[i+3].split()[3]), float(self.lines[i+3].split()[4]), float(self.lines[i+3].split()[5])])
            elif self.lines[i].split()[0] == "site" and self.lines[i].split()[-1] == "units)" and self.lines[i].split()[-2] == "(alat":
                self.run_info["site_line_number"] = i
            elif self.lines[i].split()[0] == "number" and self.lines[i].split()[2] == 'k':
                if self.lines[i].split()[5] == "(tetrahedron":
                    self.scf_params["degauss"] = "tetrahedron method: degauss not needed"
                else:
                    self.scf_params["degauss"] = float(self.lines[i].split()[9])
                self.run_info["number-of-k-points"] = int(self.lines[i].split()[4])
            elif self.lines[i].split()[0] == "Estimated" and self.lines[i].split()[1] == "max":
                self.run_info["ram_per_process"] = self.lines[i].split()[7] + " " + self.lines[i].split()[8]
                self.run_info["total_ram"] = self.lines[i+2].split()[5] + " " + self.lines[i+2].split()[6]
            elif self.lines[i].split()[0] == "total" and self.lines[i].split()[1] == "energy":
                # the total energy of the last iteration is not print like the previous scf iteration
                # it begin with a ! total energy
                self.run_info["scf_energies"].append(float(self.lines[i].split()[3]))
            elif self.lines[i].split()[0] == "!" and self.lines[i].split()[5] == "Ry":
                self.run_info["scf_final_energy"] = float(self.lines[i].split()[4])
            elif self.lines[i].split()[0] == "convergence" and self.lines[i].split()[3] == "achieved":
                self.run_info["scf_iterations"] = int(self.lines[i].split()[5])
            elif self.lines[i].split()[0] ==  "the" and self.lines[i].split()[1] == "Fermi":
                self.run_info["fermi_energy"] = float(self.lines[i].split()[4])
            elif self.lines[i].split()[0] == "Total" and self.lines[i].split()[1] == "force":
                self.run_info["total_force"] = float(self.lines[i].split()[3])
            elif self.lines[i].split()[0] == "Computing" and self.lines[i].split()[-1] == "pressure":
                self.run_info["total_stress_ry_bohr_3"] = []
                self.run_info["total_stress_kbar"] = []
                self.run_info["pressure"] = float(self.lines[i+2].split()[-1])
                self.run_info["total_stress_ry_bohr_3"].append([float(self.lines[i+3].split()[0]), float(self.lines[i+3].split()[1]), float(self.lines[i+3].split()[2])])
                self.run_info["total_stress_ry_bohr_3"].append([float(self.lines[i+4].split()[0]), float(self.lines[i+4].split()[1]), float(self.lines[i+4].split()[2])])
                self.run_info["total_stress_ry_bohr_3"].append([float(self.lines[i+5].split()[0]), float(self.lines[i+5].split()[1]), float(self.lines[i+5].split()[2])])
                self.run_info["total_stress_kbar"].append([float(self.lines[i+3].split()[3]), float(self.lines[i+3].split()[4]), float(self.lines[i+3].split()[5])])
                self.run_info["total_stress_kbar"].append([float(self.lines[i+4].split()[3]), float(self.lines[i+4].split()[4]), float(self.lines[i+4].split()[5])])
                self.run_info["total_stress_kbar"].append([float(self.lines[i+5].split()[3]), float(self.lines[i+5].split()[4]), float(self.lines[i+5].split()[5])])

        # note: at present len(self.run_info["scf_energies"]) = len(self.run_info["scf_iterations"]) - 1
        # because the total energy of the last step is not printed in format like the previous scf step,
        # and it is printed as the '!    total energy              = ' where there is a "!" in the beginning
        # now we append the final scf step energy to self.run_info["scf_energies"]
        self.run_info["scf_energies"].append(self.run_info["scf_final_energy"])
        # ----------------------------------------------------------------------
        # get the xyz structure from information extracted above:
        self.xyz = base_xyz()
        self.xyz.natom = self.scf_params["nat"]
        begin = self.run_info["site_line_number"] + 1
        # Warning:
        # there are numeric erros when obtaining atom coordinated from qe output
        # in unit of alat and multiplied by alat and bohr. namely the acquired
        # atomic coordinates have numeric errors compared to the input xyz
        # so be cautious when use it.
        bohr = 0.529177208   # 1 Bohr = 0.529177208 Angstrom
        for i in range(self.xyz.natom):
            self.xyz.atoms.append(Atom(
                self.lines[begin+i].split()[1],
                self.scf_params["alat_au"] * bohr * float(self.lines[begin+i].split()[6]),
                self.scf_params["alat_au"] * bohr * float(self.lines[begin+i].split()[7]),
                self.scf_params["alat_au"] * bohr * float(self.lines[begin+i].split()[8])))
        self.xyz.cell = self.scf_params["cell_a_alat"] # now in unit of alat

        for i in range(3):
            for j in range(3):
                self.xyz.cell[i][j] = self.scf_params["cell_a_alat"][i][i] * self.scf_params["alat_au"] * bohr
        # now self.xyz.cell are in unit of Angstrom
        # ----------------------------------------------------------------------

    # def

# stronger please
