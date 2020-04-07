#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import copy
import json
import datetime
import subprocess
import matplotlib.pyplot as plt

from pymatflow.base.xyz import base_xyz
from pymatflow.base.atom import Atom


class opt_out:
    """
    """
    def __init__(self):
        """
        output is the output file of scf run
        """
        self.file = None
        self.opt_params = {}
        self.run_info = {}

        self.run_type = None
        self.job_done = None # whether calculation has finished
        self.relaxed = None # whether structure is relaxed or vc-relaxed

        self.cell = None # optimized cell now only used when run_type == vc-relax
        self.trajectory = None

    def get_info(self, file):
        """
        get the general information from relax or vc-relax run output file
        which is now stored in self.lines
        """
        self.clean()

        self.file = file
        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()

        # check whether calculation is finished
        if len(self.lines[-2].split()) == 2 and self.lines[-2].split()[0] == "JOB" and self.lines[-2].split()[1] == "DONE.":
            self.job_done = True
        else:
            self.job_done = False

        # check the run_type: relax or vc-relax
        self.run_type = 'relax'
        for line in self.lines:
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "CELL_PARAMETERS" and line.split()[1].split("\n")[0] == "(angstrom)":
                self.run_type = 'vc-relax'
                break

        # check whether successfully relaxed
        self.relaxed = False
        for line in self.lines:
            if line == "Begin final coordinates\n":
                self.relaxed = True
                break

        self.get_trajectory()
        self.get_opt_params_and_run_info()

    def clean(self):
        self.file = None
        self.opt_params = {}
        self.run_info = {}

    #
    def get_opt_params_and_run_info(self):
        """
        self.run_info[]
            start_time: the task start time
            stop_time: the task stop time
            #scf_energies: all the energies during the scf procedure
            #fermi_energy: fermi energy of the system (if output)

        """
        #self.run_info["scf_energies"] = []

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
                self.opt_params["alat_au"] = float(self.lines[i+1].split()[4])
                self.opt_params["nat"] = int(self.lines[i+3].split()[4])
                self.opt_params["nelectron"] = float(self.lines[i+5].split()[4])
                self.opt_params["n_ks_state"] = int(self.lines[i+6].split("=")[1])
                self.opt_params["ecutwfc"] = int(float(self.lines[i+7].split()[3]))
                self.opt_params["ecutrho"] = int(float(self.lines[i+8].split()[4]))
                self.opt_params["conv_thr"] = float(self.lines[i+9].split()[3])
                self.opt_params["mixing_beta"] = float(self.lines[i+10].split()[3])
                if "nstep" in self.opt_params and self.opt_params["nstep"] != None:
                    pass
                else:
                    self.opt_params["nstep"] = int(self.lines[i+13].split()[2])
            elif self.lines[i].split()[0] == "crystal" and self.lines[i].split()[1] == "axes:" and self.lines[i].split()[-1] =="alat)":
                self.opt_params["cell_a_alat"] = []
                self.opt_params["cell_a_alat"].append([float(self.lines[i+1].split()[3]), float(self.lines[i+1].split()[4]), float(self.lines[i+1].split()[5])])
                self.opt_params["cell_a_alat"].append([float(self.lines[i+2].split()[3]), float(self.lines[i+2].split()[4]), float(self.lines[i+2].split()[5])])
                self.opt_params["cell_a_alat"].append([float(self.lines[i+3].split()[3]), float(self.lines[i+3].split()[4]), float(self.lines[i+3].split()[5])])
            elif self.lines[i].split()[0] == "reciprocal" and self.lines[i].split()[1] == "axes:" and self.lines[i].split()[-1] == "pi/alat)": # actually '2 pi/alat'
                self.opt_params["cell_b_2pi_alat"] = []
                self.opt_params["cell_b_2pi_alat"].append([float(self.lines[i+1].split()[3]), float(self.lines[i+1].split()[4]), float(self.lines[i+1].split()[5])])
                self.opt_params["cell_b_2pi_alat"].append([float(self.lines[i+2].split()[3]), float(self.lines[i+2].split()[4]), float(self.lines[i+2].split()[5])])
                self.opt_params["cell_b_2pi_alat"].append([float(self.lines[i+3].split()[3]), float(self.lines[i+3].split()[4]), float(self.lines[i+3].split()[5])])
            elif self.lines[i].split()[0] == "site" and self.lines[i].split()[-1] == "units)" and self.lines[i].split()[-2] == "(alat":
                self.run_info["site_line_number"] = i
            elif self.lines[i].split()[0] == "number" and self.lines[i].split()[2] == 'k':
                if self.lines[i].split()[5] == "(tetrahedron":
                    self.opt_params["degauss"] = "tetrahedron method: degauss not needed"
                else:
                    self.opt_params["degauss"] = float(self.lines[i].split()[9])
                self.run_info["number-of-k-points"] = int(self.lines[i].split()[4])
            elif self.lines[i].split()[0] == "Estimated" and self.lines[i].split()[1] == "max":
                self.run_info["ram_per_process"] = self.lines[i].split()[7] + " " + self.lines[i].split()[8]
                self.run_info["total_ram"] = self.lines[i+2].split()[5] + " " + self.lines[i+2].split()[6]


        # ----------------------------------------------------------------------
        # get the input xyz structure from information extracted above:
        self.xyz = base_xyz()
        self.xyz.natom = self.opt_params["nat"]
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
                self.opt_params["alat_au"] * bohr * float(self.lines[begin+i].split()[6]),
                self.opt_params["alat_au"] * bohr * float(self.lines[begin+i].split()[7]),
                self.opt_params["alat_au"] * bohr * float(self.lines[begin+i].split()[8])))
        self.xyz.cell = self.opt_params["cell_a_alat"] # now in unit of alat

        for i in range(3):
            for j in range(3):
                self.xyz.cell[i][j] = self.opt_params["cell_a_alat"][i][i] * self.opt_params["alat_au"] * bohr
        # now self.xyz.cell are in unit of Angstrom
        # ----------------------------------------------------------------------

        # get information output each ion step
        self._get_info_for_each_ions_step()

    def _get_info_for_each_ions_step(self):
        self.run_info["total_energy_each_ion_step"] = []
        self.run_info["fermi_energy_each_ion_step"] = []
        self.run_info["total_force_each_ion_step"] = []
        self.run_info["total_force_scf_correction_each_ion_step"] = []
        self.run_info["scf_iterations_each_ion_step"] = []
        for i in range(len(self.lines)):
            # if it is an empty line continue to next line
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0] == "!" and self.lines[i].split()[1] == "total" and self.lines[i].split()[2] == "energy":
                self.run_info["total_energy_each_ion_step"].append(float(self.lines[i].split()[4])) # in unit of Ry
            if self.lines[i].split()[0] == "the" and self.lines[i].split()[1] == "Fermi":
                self.run_info["fermi_energy_each_ion_step"].append(float(self.lines[i].split()[4])) # in unit of eV
            if self.lines[i].split()[0] == "Total" and self.lines[i].split()[1] == "force" and self.lines[i].split()[2] == "=":
                self.run_info["total_force_each_ion_step"].append(float(self.lines[i].split()[3])) # in unit of Ry/au
                self.run_info["total_force_scf_correction_each_ion_step"].append(self.lines[i].split()[8]) # in unit of Ry/au
            if self.lines[i] .split()[0] == "convergence" and self.lines[i].split()[1] == "has":
                self.run_info["scf_iterations_each_ion_step"].append(int(self.lines[i].split()[5]))
            

    def get_trajectory(self):
        """
        Note: initial input structure is not in self.trajectory, but is in self.xyz
            self.trajectory contains all other structures and the last of it
            is the optimized structure.
        """
        self.trajectory = []
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "ATOMIC_POSITIONS":
                atm = []
                j = i + 1
                while len(self.lines[j].split()) == 4:
                    atm.append(Atom(self.lines[j].split()[0], float(self.lines[j].split()[1]), float(self.lines[j].split()[2]), float(self.lines[j].split()[3])))
                    j = j + 1
                self.trajectory.append(atm)
        #
        if self.relaxed == True and self.run_type == "vc-relax":
            # get the line number of the 'Begin final coordinates'
            # and 'End final coordinates'
            begin_final_coord_line = 0
            end_final_coord_line = 0
            while self.lines[begin_final_coord_line] != "Begin final coordinates\n":
                begin_final_coord_line += 1
            while self.lines[end_final_coord_line] != "End final coordinates\n":
                end_final_coord_line += 1
            # get the optimized cell
            self.final_cell = []
            for i in range(begin_final_coord_line+5, begin_final_coord_line+8):
                vec = []
                for j in range(3):
                    vec.append(float(self.lines[i].split()[j]))
                self.final_cell.append(vec)
        #
    def print_final_structure(self, xyz="optimized.xyz"):
        if self.relaxed == False:
            with open("final-structure(not-relaxed).xyz", 'w') as fout:
                fout.write("%d\n" % len(self.trajectory[0]))
                fout.write("Warning(%s): structure failed to be relaxed or vc-relaxed, this is the final structure(unrelaxed)\n" % self.run_type)
                for atom in self.trajectory[-1]:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
            return
        # printout relaxed structure
        with open(xyz, 'w') as fout:
            fout.write("%d\n" % len(self.trajectory[-1]))
            if self.run_type == "vc-relax":
                cell = self.final_cell
                fout.write("cell: %.9f %.9f %.9f | %.9f %.9f %.9f | %.9f %.9f %.9f\n" % (cell[0][0], cell[0][1], cell[0][2], cell[1][0], cell[1][1], cell[1][2], cell[2][0], cell[2][1], cell[2][2]))
            else:
                cell = self.xyz.cell
                fout.write("cell: %.9f %.9f %.9f | %.9f %.9f %.9f | %.9f %.9f %.9f\n" % (cell[0][0], cell[0][1], cell[0][2], cell[1][0], cell[1][1], cell[1][2], cell[2][0], cell[2][1], cell[2][2]))
            for atom in self.trajectory[-1]:
                fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))

    def print_trajectory(self, xyz="trajectory.xyz"):
        with open(xyz, 'w') as fout:
            for i in range(len(self.trajectory)):
                fout.write("%d\n" % len(self.trajectory[i]))
                fout.write("i = %d\n" % i)
                for atom in self.trajectory[i]:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))


    def plot_run_info(self):
        """
        """
        plt.plot(self.run_info["scf_iterations_each_ion_step"])
        plt.title("Iterations each ion step")
        plt.xlabel("Ion step")
        plt.ylabel("Scf iterations")
        plt.tight_layout()
        plt.savefig("scf-iterations-each-ion-step.png")
        plt.close()

        plt.plot(self.run_info["total_energy_each_ion_step"])
        plt.title("Total energies each ion step")
        plt.xlabel("Ion step")
        plt.ylabel("Total Energy (Ry)")
        plt.tight_layout()
        plt.savefig("total-energy-each-ion-step.png")
        plt.close()

        plt.plot(self.run_info["fermi_energy_each_ion_step"])
        plt.title("Fermi energies each ion step")
        plt.xlabel("Ion step")
        plt.ylabel("Fermi energiy (eV)")
        plt.tight_layout()
        plt.savefig("fermi-energy-each-ion-step.png")
        plt.close()

        plt.plot(self.run_info["total_force_each_ion_step"])
        plt.title("Total force each ion step")
        plt.xlabel("Ion step")
        plt.ylabel("Total forces (Ry/au)")
        plt.tight_layout()
        plt.savefig("total-force-each-ion-step.png")
        plt.close()

    def markdown_report(self, md="opt-info.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# 几何优化实验统计\n")
            fout.write("几何优化类型: %s\n" % self.run_type)
            fout.write("几何优化任务是否结束:%s\n" % str(self.job_done))
            if self.job_done == True:
                fout.write("是否成功优化: %s\n" % str(self.relaxed))
            else:
                fout.write("是否成功优化: %s\n" % ("运行未结束, 结果未知"))
            fout.write("## 优化参数\n")
            for item in self.opt_params:
                fout.write("- %s: %s\n" % (item, str(self.opt_params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out
            # Importante: the length of the time string might be different, depending
            # on the value of hours and minutes and seconds. if they are two digits
            # number, they will be divided like: '11: 6: 2', only when they all are
            # two digtis number, they will not be divided '11:16:12'
            # so we have to preprocess it to build the right time string to pass into
            # datetime.datetime.strptime()
            if len(self.run_info["start_time"].split()) == 8:
                start_str = self.run_info["start_time"].split()[5]+"-"+self.run_info["start_time"].split()[7]
            elif len(self.run_info["start_time"].split()) == 9:
                start_str = self.run_info["start_time"].split()[5]+"-"+self.run_info["start_time"].split()[7]+self.run_info["start_time"].split()[8]
            elif len(self.run_info["start_time"].split()) == 10:
                start_str = self.run_info["start_time"].split()[5]+"-"+self.run_info["start_time"].split()[7]+self.run_info["start_time"].split()[8]+self.run_info["start_time"].split()[9]
            else:
                print("===============================================\n")
                print("                  Warning !!!\n")
                print("===============================================\n")
                print("qe.post.opt.markdown_report:\n")
                print("failed to parse start_time string\n")
                sys.exit(1)
            if self.job_done == True:
                if len(self.run_info["stop_time"].split()) == 7:
                    stop_str = self.run_info["stop_time"].split()[6]+"-"+self.run_info["stop_time"].split()[5]
                elif len(self.run_info["stop_time"].split()) == 8:
                    stop_str = self.run_info["stop_time"].split()[7]+"-"+self.run_info["stop_time"].split()[5]+self.run_info["stop_time"].split()[6]
                elif len(self.run_info["stop_time"].split()) == 9:
                    stop_str = self.run_info["stop_time"].split()[8]+"-"+self.run_info["stop_time"].split()[5]+self.run_info["stop_time"].split()[6]+self.run_info["stop_time"].split()[7]
                else:
                    print("===============================================\n")
                    print("                  Warning !!!\n")
                    print("===============================================\n")
                    print("qe.post.opt.markdown_report:\n")
                    print("failed to parse stop_time string\n")
                    sys.exit(1)

            start = datetime.datetime.strptime(start_str, "%d%b%Y-%H:%M:%S")
            if self.job_done == True:
                stop = datetime.datetime.strptime(stop_str, "%d%b%Y-%H:%M:%S")
                delta_t = stop -start
            fout.write("- Time consuming:\n")
            if self.job_done == True:
                fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
            else:
                fout.write("  - job is not finished yet, but it starts at %s\n" % start)
            # end the time information
            for item in self.run_info:
                fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

            fout.write("## 运行信息图示\n")
            fout.write("Iterations per SCF\n")
            fout.write("![Iterations per SCF](scf-iterations-each-ion-step.png)\n")

            fout.write("Total energies per SCF\n")
            fout.write("![Total energies per SCF](total-energy-each-ion-step.png)\n")

            fout.write("Fermi energies per SCF\n")
            fout.write("![Fermi energies per SCF](fermi-energy-each-ion-step.png)\n")

            fout.write("Total forces per SCF\n")
            fout.write("![Total forces per SCF](total-force-each-ion-step.png)\n")

    def to_json_string(self):
        """
        Note:
            self.run_info is a dict, but it cannot be dumped to json directory by json.dumps
            because there are objects inside self.run_info, like datetime and Atom() that cannot
            be serialized by json.dumps()
        :return out: json string processed by json.dumps()
        """
        run_info = copy.deepcopy(self.run_info)
        opt_params = copy.deepcopy(self.opt_params)
        # convert datetime to str
        run_info["start_time"] = str(run_info["start_time"])
        # merge opt_params and run_info and output to json and return
        out = opt_params
        out.update(run_info)
        return json.dumps(out, indent=4)


    def export(self, directory):
        """
        Note:
            * will only printout the final structure if the job is done
        """
        os.chdir(directory)
        os.system("mkdir -p post-processing")
        os.chdir("post-processing")
        if self.job_done == True:
            self.print_final_structure()
        self.print_trajectory()
        self.plot_run_info()
        self.markdown_report("opt-info.md")
        with open("opt.json", 'w') as fout:
            fout.write(self.to_json_string())
        os.system("cd ../../")


class opt_post:
    """
    Note:
        opt_post can extract information for the geometric optimization running,
        including 'relax' and 'vc-relax'.
        it will printout the trajectory file(xyz format) and the final optimized
        structure(if not relaxed, the final structure of the running).
        so even when your ion step not converged within maximum steps, it can also
        extract the structure in the final ion step and print it out.

        plus now opt_post can also process the output file of geometric optimziation
        even if the job is not finished yet and generate the Running information.

        we must know when it is a relax geometric optimization, the number of scf
        cycles equals to the length of the trajectory. but when it is a vc-relax
        running, the program will do a final scf calculation on the relaxed structure
        , and during that calculation it will not print out structure information.
        so the length of trajectory of a vc-relax running might be one less than the
        number of scf cycles.(will be euqal if the vc-relax is unrelxed, so there won't
        be a final scf)
    """
    def __init__(self, output, run_type):
        """
        output:
            the output file of opt run. it could be a geo opt or a cell opt

        relaxed:
            whether the structure successfully relaxed.
        """
        self.file = output
        self.run_type = run_type
        self.opt_params = {}
        self.run_info = {}
        self.job_done = None # whether calculation has finished
        self.relaxed = None # whether structure is relaxed or vc-relaxed

        self.cell = None # optimized cell now only used when run_type == vc-relax
        self.trajectory = None

        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()
        self.get_info()

    def get_info(self):
        """
        get the general information of opt run from opt run output file
        which is now stored in self.lines
        """
        # check whether calculation is finished
        if len(self.lines[-2].split()) == 2 and self.lines[-2].split()[0] == "JOB" and self.lines[-2].split()[1] == "DONE.":
            self.job_done = True
        else:
            self.job_done = False
        # check whether successfully relaxed
        self.relaxed = False
        for line in self.lines:
            if line == "Begin final coordinates\n":
                self.relaxed = True
                break

        self.get_trajectory()
        self.get_opt_params_and_run_info()


    def get_trajectory(self):
        self.trajectory = []
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "ATOMIC_POSITIONS":
                atm = []
                j = i + 1
                while len(self.lines[j].split()) == 4:
                    atm.append(Atom(self.lines[j].split()[0], float(self.lines[j].split()[1]), float(self.lines[j].split()[2]), float(self.lines[j].split()[3])))
                    j = j + 1
                self.trajectory.append(atm)
        #
        if self.relaxed == True and self.run_type == "vc-relax":
            # get the line number of the 'Begin final coordinates'
            # and 'End final coordinates'
            begin_final_coord_line = 0
            end_final_coord_line = 0
            while self.lines[begin_final_coord_line] != "Begin final coordinates\n":
                begin_final_coord_line += 1
            while self.lines[end_final_coord_line] != "End final coordinates\n":
                end_final_coord_line += 1
            # get the optimized cell
            self.cell = []
            for i in range(begin_final_coord_line+5, begin_final_coord_line+8):
                for j in range(3):
                    self.cell.append(float(self.lines[i].split()[j]))
        #

    def get_opt_params_and_run_info(self):
        """
        run_info["iterations"]: scf iterations per scf step
        run_info["total-energies"]: total energies of every scf step
        run_info["fermi-energies"]: fermi energies of every scf step
        run_info["total-forces"]: total forces of every scf step
        """
        self.run_info["iterations"] = []
        self.run_info["total-energies"] = []
        self.run_info["fermi-energies"] = []
        self.run_info["total-forces"] = []

        for line in self.lines:
            # if it is an empty line continue to next line
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "Program" and line.split()[1] == "PWSCF" and line.split()[3] == "starts":
                self.run_info["start-time"] = line.split("\n")[0]
            if line.split()[0] == "This" and line.split()[1] == "run" and line.split()[3] == "terminated":
                self.run_info["stop-time"] = line.split("\n")[0]
            if line.split()[0] == "kinetic-energy":
                self.opt_params["ecutwfc"] = int(float(line.split()[3]))
            if line.split()[0] == "convergence threshold":
                self.opt_params["conv_thr"] = float(line.split()[3])
            if line.split()[0] == "mixing" and line.split()[1] == 'beta':
                self.opt_params["mixing_beta"] = float(line.split()[3])
            if line.split()[0] == "number" and line.split()[2] == 'k':
                self.opt_params["degauss"] = float(line.split()[9])
            if line.split()[0] == "convergence" and line.split()[3] == "achieved":
                self.run_info["iterations"].append(int(line.split()[5]))
            if line.split()[0] == "!" and line.split()[5] == "Ry":
                self.run_info["total-energies"].append(float(line.split()[4]))
            if line.split()[0] ==  "the" and line.split()[1] == "Fermi":
                self.run_info["fermi-energies"].append(float(line.split()[4]))
            if line.split()[0] == "Total" and line.split()[1] == "force":
                self.run_info["total-forces"].append(float(line.split()[3]))

        self.run_info["scf-cycles"] = len(self.run_info["iterations"])
        if self.run_type == "relax":
            self.run_info["ion-steps"] = len(self.run_info["iterations"]) - 1
        elif self.run_type == "vc-relax":
            self.run_info["ion-steps"] = len(self.run_info["iterations"]) - 2

    def print_final_structure(self, xyz="optimized.xyz"):
        if self.relaxed == False:
            with open("final-structure(not-relaxed).xyz", 'w') as fout:
                fout.write("%d\n" % len(self.trajectory[0]))
                fout.write("Warning(%s): structure failed to be relaxed or vc-relaxed, this is the final structure(unrelaxed)\n" % self.run_type)
                for atom in self.trajectory[-1]:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
            return
        # printout relaxed structure
        cell = self.cell
        with open(xyz, 'w') as fout:
            fout.write("%d\n" % len(self.trajectory[-1]))
            if self.run_type == "vc-relax":
                fout.write("cell: %.9f %.9f %.9f | %.9f %.9f %.9f | %.9f %.9f %.9f\n" % (cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], cell[6], cell[7], cell[8]))
            else:
                fout.write("type of opt run: relax -> the cell is not changed, so go and find the original cell\n")
            for atom in self.trajectory[-1]:
                fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))

    def print_trajectory(self, xyz="trajectory.xyz"):
        with open(xyz, 'w') as fout:
            for i in range(len(self.trajectory)):
                fout.write("%d\n" % len(self.trajectory[i]))
                fout.write("i = %d\n" % i)
                for atom in self.trajectory[i]:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))

    def view_trajectory(self, trajfile="trajectory.xyz"):
        #os.system("xcrysden --xyz %s" % trajfile)
        subprocess.call(["xcrysden", "--xyz", trajfile])

    def plot_run_info(self):
        """
        """
        plt.plot(self.run_info["iterations"])
        plt.title("Iterations per SCF")
        plt.xlabel("Scf cycles")
        plt.ylabel("iterations")
        plt.tight_layout()
        plt.savefig("iterations-per-scf.png")
        plt.close()

        plt.plot(self.run_info["total-energies"])
        plt.title("Total energies per SCF")
        plt.xlabel("Scf cycles")
        plt.ylabel("Total Energies (Ry)")
        plt.tight_layout()
        plt.savefig("total-energies-per-scf.png")
        plt.close()

        plt.plot(self.run_info["fermi-energies"])
        plt.title("Fermi energies per SCF")
        plt.xlabel("Scf cycles")
        plt.ylabel("Fermi energies (eV)")
        plt.tight_layout()
        plt.savefig("fermi-energies-per-scf.png")
        plt.close()

        plt.plot(self.run_info["total-forces"])
        plt.title("Total forces per SCF")
        plt.xlabel("Scf cycles")
        plt.ylabel("Total forces (Ry/au)")
        plt.tight_layout()
        plt.savefig("total-forces-per-scf.png")
        plt.close()

    def markdown_report(self, md="OptimizationReport.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# 几何优化实验统计\n")
            fout.write("几何优化类型: %s\n" % self.run_type)
            fout.write("几何优化任务是否结束:%s\n" % str(self.job_done))
            if self.job_done == True:
                fout.write("是否成功优化: %s\n" % str(self.relaxed))
            else:
                fout.write("是否成功优化: %s\n" % ("运行未结束, 结果未知"))
            fout.write("## 优化参数\n")
            for item in self.opt_params:
                fout.write("- %s: %s\n" % (item, str(self.opt_params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out
            # Importante: the length of the time string might be different, depending
            # on the value of hours and minutes and seconds. if they are two digits
            # number, they will be divided like: '11: 6: 2', only when they all are
            # two digtis number, they will not be divided '11:16:12'
            # so we have to preprocess it to build the right time string to pass into
            # datetime.datetime.strptime()
            if len(self.run_info["start-time"].split()) == 8:
                start_str = self.run_info["start-time"].split()[5]+"-"+self.run_info["start-time"].split()[7]
            elif len(self.run_info["start-time"].split()) == 9:
                start_str = self.run_info["start-time"].split()[5]+"-"+self.run_info["start-time"].split()[7]+self.run_info["start-time"].split()[8]
            elif len(self.run_info["start-time"].split()) == 10:
                start_str = self.run_info["start-time"].split()[5]+"-"+self.run_info["start-time"].split()[7]+self.run_info["start-time"].split()[8]+self.run_info["start-time"].split()[9]
            else:
                print("===============================================\n")
                print("                  Warning !!!\n")
                print("===============================================\n")
                print("qe.post.opt.markdown_report:\n")
                print("failed to parse start-time string\n")
                sys.exit(1)
            if self.job_done == True:
                if len(self.run_info["stop-time"].split()) == 7:
                    stop_str = self.run_info["stop-time"].split()[6]+"-"+self.run_info["stop-time"].split()[5]
                elif len(self.run_info["stop-time"].split()) == 8:
                    stop_str = self.run_info["stop-time"].split()[7]+"-"+self.run_info["stop-time"].split()[5]+self.run_info["stop-time"].split()[6]
                elif len(self.run_info["stop-time"].split()) == 9:
                    stop_str = self.run_info["stop-time"].split()[8]+"-"+self.run_info["stop-time"].split()[5]+self.run_info["stop-time"].split()[6]+self.run_info["stop-time"].split()[7]
                else:
                    print("===============================================\n")
                    print("                  Warning !!!\n")
                    print("===============================================\n")
                    print("qe.post.opt.markdown_report:\n")
                    print("failed to parse stop-time string\n")
                    sys.exit(1)

            start = datetime.datetime.strptime(start_str, "%d%b%Y-%H:%M:%S")
            if self.job_done == True:
                stop = datetime.datetime.strptime(stop_str, "%d%b%Y-%H:%M:%S")
                delta_t = stop -start
            fout.write("- Time consuming:\n")
            if self.job_done == True:
                fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
            else:
                fout.write("  - job is not finished yet, but it starts at %s\n" % start)
            # end the time information
            for item in self.run_info:
                fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

            fout.write("## 运行信息图示\n")
            fout.write("Iterations per SCF\n")
            fout.write("![Iterations per SCF](iterations-per-scf.png)\n")

            fout.write("Total energies per SCF\n")
            fout.write("![Total energies per SCF](total-energies-per-scf.png)\n")

            fout.write("Fermi energies per SCF\n")
            fout.write("![Fermi energies per SCF](fermi-energies-per-scf.png)\n")

            fout.write("Total forces per SCF\n")
            fout.write("![Total forces per SCF](total-forces-per-scf.png)\n")


    def export(self):
        """
        Note:
            * will only printout the final structure if the job is done
        """
        if self.job_done == True:
            self.print_final_structure()
        self.print_trajectory()
        self.plot_run_info()
        self.markdown_report("OptimizationReport.md")
