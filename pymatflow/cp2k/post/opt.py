#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import copy
import json
import datetime
import subprocess
import matplotlib.pyplot as plt

from pymatflow.cmd.structflow import read_structure, write_structure


class OptOut:
    """
    """
    def __init__(self):
        """
        output is the output file of scf run
        """
        self.file = None
        self.params = {}
        self.run_info = {}
        self.run_type = None # GEO_OPT / CELL_OPT

    def get_info(self, file):
        """
        get the general information of running from opt run output file
        which is now stored in self.lines
        """
        self.clean()

        self.file = file
        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()
        self.get_params_and_run_info()


    def clean(self):
        self.file = None
        self.params = {}
        self.run_info = {}

    #
    def get_params_and_run_info(self):
        """
        self.run_info[]
            start_time: the task start time
            stop_time: the task stop time
            scf_energies: all the energies during the scf procedure
            #fermi_energy: fermi energy of the system (if output)

        """
        self.run_info["total_energy_each_ion_step"] = []

        for i in range(len(self.lines)):
            # if it is an empty line continue to next line
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0] == "****" and self.lines[i].split()[1] == "****" and self.lines[i].split()[5] == "STARTED":
                self.run_info["start_time"] = self.lines[i].split("\n")[0]
            elif self.lines[i].split()[0] == "****" and self.lines[i].split()[1] == "****" and self.lines[i].split()[5] == "ENDED":
                self.run_info["stop_time"] = self.lines[i].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Basis":
                self.params["BASIS_SET_FILE_NAME"] = self.lines[i].split()[-1].split("\n")[0]
                self.params["POTENTIAL_FILE_NAME"] = self.lines[i+1].split()[-1].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Project":
                self.params["PROJECT_NAME"] = self.lines[i].split()[-1].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Run":
                self.params["RUN_TYPE"] = self.lines[i].split()[-1].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Global":
                self.params["PRINT_LEVEL"] = self.lines[i].split()[-1].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Total" and self.lines[i].split()[2] == "number":
                self.run_info["mpi_processes"] = int(self.lines[i].split()[-1].split("\n")[0])
                self.run_info["cpu_model"] = self.lines[i+3].split("name")[1].split("\n")[0]
            elif self.lines[i].split()[0] == "-" and self.lines[i].split()[1] == "Atoms:":
                self.run_info["n_atom"] = int(self.lines[i].split()[-1].split("\n")[0])
                self.run_info["n_shell"] = int(self.lines[i+2].split()[-1].split("\n")[0])
            elif self.lines[i].split()[0] == "max_scf:":
                self.params["MAX_SCF"] = int(self.lines[i].split()[-1].split("\n")[0])
                self.params["MAX_SCF_HISTORY"] = int(self.lines[i+1].split()[-1].split("\n")[0])
                self.params["MAX_DIIS"] = int(self.lines[i+2].split()[-1].split("\n")[0])
            elif self.lines[i].split()[0] == "eps_scf:":
                self.params["EPS_SCF"] = float(self.lines[i].split()[1])
                self.params["EPS_SCF_HISTORY"] = float(self.lines[i+1].split()[1])
                self.params["EPS_DIIS"] = float(self.lines[i+2].split()[1])
            elif self.lines[i].split()[0] == "Mixing" and self.lines[i].split()[1] == 'method:':
                self.params["MIXING"] = self.lines[i].split()[2]
            elif self.lines[i].split()[0] == "added" and self.lines[i].split()[1] == "MOs":
                self.params["ADDED_MOS"] = int(self.lines[i].split()[2])
            elif self.lines[i].split()[0] == "Number" and self.lines[i].split()[1] == "of" and self.lines[i].split()[2] == "electrons:":
                self.run_info["n_electrons"] = int(self.lines[i].split()[3])
                self.run_info["n_occcupied_orbital"] = int(self.lines[i+1].split()[4])
                self.run_info["n_molecular_orbital"] = int(self.lines[i+2].split()[4])
                self.run_info["n_orbital_function"] = int(self.lines[i+4].split()[4])
            elif self.lines[i].split()[0] == "ENERGY|" and self.lines[i].split()[4] == "QS":
                self.run_info["total_energy_each_ion_step"].append(float(self.lines[i].split()[8])) # in unit of a.u.
            else:
                pass
        # ----------------------------------------------------------------------
        # get the xyz structure from information extracted above:
        # WARNING: in low level print of cp2k, there is no structure coordinates
        # in the output file
        # ----------------------------------------------------------------------
        #
    def plot_info(self):
        # now output the scf information in the current directory(post-processing)
        plt.plot(self.run_info["total_energy_each_ion_step"], marker="o")
        plt.title("Total energy each ion step")
        plt.xlabel("Ion step")
        plt.ylabel("Total energy(a.u.)")
        plt.tight_layout()
        plt.savefig("total-energy-each-ion-step.png")
        plt.close()

    def get_optimized_structure_geo_opt(self, xyztraj, geo_opt_inp, directory):
        """
        :param xyztraj: output xxx.xyz trajectory file
        :param geo_opt_inp: input file for geo opt
        :param directory: directory to put the optimized structure

        Note: deal with GEO_OPT
        """

        with open(geo_opt_inp, 'r') as fin:
            geo_opt_inp_lines = fin.readlines()
        for i in range(len(geo_opt_inp_lines)):
            if len(geo_opt_inp_lines[i].split()) > 0 and geo_opt_inp_lines[i].split()[0].upper() == "&CELL":
                j = 1
                while not ("&END" in geo_opt_inp_lines[i+j].upper().split() and "CELL" in geo_opt_inp_lines[i+j].upper().split()):
                    if len(geo_opt_inp_lines[i+j].split()) > 0 and geo_opt_inp_lines[i+j].split()[0].upper() == "A":
                        a = []
                        for k in range(3):
                            a.append(float(geo_opt_inp_lines[i+j].split()[k+1]))
                    elif len(geo_opt_inp_lines[i+j].split()) > 0 and geo_opt_inp_lines[i+j].split()[0].upper() == "B":
                        b = []
                        for k in range(3):
                            b.append(float(geo_opt_inp_lines[i+j].split()[k+1]))
                    elif len(geo_opt_inp_lines[i+j].split()) > 0 and geo_opt_inp_lines[i+j].split()[0].upper() == "C":
                        c = []
                        for k in range(3):
                            c.append(float(geo_opt_inp_lines[i+j].split()[k+1]))
                    j += 1
                break
        cell = []
        cell.append(a)
        cell.append(b)
        cell.append(c)
        
        traj_lines = []
        with open(xyztraj, 'r') as fin:
            for line in fin:
                if len(line.split()) == 0:
                    continue
                traj_lines.append(line)
        natom = int(traj_lines[0].split()[0])
        nimage = int(len(traj_lines)/(natom+2))
        with open(os.path.join(directory, "geo-opt-optimized.xyz"), 'w') as fout:
            fout.write("%d\n" %natom)
            fout.write("cell: %.9f %.9f %.9f | %.9f %.9f %.9f | %.9f %.9f %.9f\n" % (
                cell[0][0], cell[0][1], cell[0][2],
                cell[1][0], cell[1][1], cell[1][2],
                cell[2][0], cell[2][1], cell[2][2],
            ))
            for i in range((natom+2)*(nimage-1)+2, (natom+2)*nimage):
                fout.write(traj_lines[i])

        a = read_structure(filepath=os.path.join(directory, "geo-opt-optimized.xyz"))
        write_structure(structure=a, filepath=os.path.join(directory, "geo-opt-optimized.cif"))
        # end get optimized structure

    def get_optimized_structure_cell_opt(self, xyztraj, cell_opt_out, directory):
        """
        :param xyztraj: output xxx.xyz trajectory file
        :param cell_opt_out: output file of cell opt
        :param directory: directory to put the optimized structure
        
        Note: deal with CELL_OPT
        """

        with open(cell_opt_out, 'r') as fin:
            cell_opt_out_lines = fin.readlines()
        for i in range(len(cell_opt_out_lines)):
            #if len(cell_opt_out_lines[i].split()) == 0:
            #    continue
            if "GEOMETRY OPTIMIZATION COMPLETED" in cell_opt_out_lines[i]:
                a = []
                b = []
                c = []
                for j in range(3):
                    a.append(float(cell_opt_out_lines[i+6].split()[4+j]))
                    b.append(float(cell_opt_out_lines[i+7].split()[4+j]))
                    c.append(float(cell_opt_out_lines[i+8].split()[4+j]))

        cell = []
        cell.append(a)
        cell.append(b)
        cell.append(c)
        
        traj_lines = []
        with open(xyztraj, 'r') as fin:
            for line in fin:
                if len(line.split()) == 0:
                    continue
                traj_lines.append(line)
        natom = int(traj_lines[0].split()[0])
        nimage = int(len(traj_lines)/(natom+2))
        with open(os.path.join(directory, "cell-opt-optimized.xyz"), 'w') as fout:
            fout.write("%d\n" %natom)
            fout.write("cell: %.9f %.9f %.9f | %.9f %.9f %.9f | %.9f %.9f %.9f\n" % (
                cell[0][0], cell[0][1], cell[0][2],
                cell[1][0], cell[1][1], cell[1][2],
                cell[2][0], cell[2][1], cell[2][2],
            ))
            for i in range((natom+2)*(nimage-1)+2, (natom+2)*nimage):
                fout.write(traj_lines[i])

        a = read_structure(filepath=os.path.join(directory, "cell-opt-optimized.xyz"))
        write_structure(structure=a, filepath=os.path.join(directory, "cell-opt-optimized.cif"))
        # end get optimized structure


    def markdown_report(self, md="opt-info.md"):
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# Optimization实验统计\n")
            fout.write("集合优化类型:%s\n" % self.params["RUN_TYPE"])
            fout.write("## 模拟参数\n")
            for item in self.params:
                fout.write("- %s: %s\n" % (item, str(self.params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out
            # note that the accuracy for the seconds is not fully guranteed
            # e.g. 2019-11-26 12:09:36.487 is read as 2019-11-26 12:09:36
            start = datetime.datetime.strptime(self.run_info["start_time"].split()[7]+"-"+self.run_info["start_time"].split()[8].split(".")[0], "%Y-%m-%d-%H:%M:%S")
            stop = datetime.datetime.strptime(self.run_info["stop_time"].split()[7]+"-"+self.run_info["stop_time"].split()[8].split(".")[0], "%Y-%m-%d-%H:%M:%S")
            delta_t = stop -start
            fout.write("- Time consuming:\n")
            fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
            # end the time information
            for item in self.run_info:
                fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

            fout.write("## 运行信息图示\n")
            fout.write("![total-energy-each-ion-step](total-energy-each-ion-step.png)")
    
    def to_json_string(self):
        """
        Note:
            self.run_info is a dict, but it cannot be dumped to json directory by json.dumps
            because there are objects inside self.run_info, like datetime that cannot
            be serialized by json.dumps()
        :return out: json string processed by json.dumps()
        """
        run_info = copy.deepcopy(self.run_info)
        params = copy.deepcopy(self.params)
        # convert datetime to str
        run_info["start_time"] = str(run_info["start_time"])
        # merge opt_params and run_info and output to json and return
        out = params
        out.update(run_info)
        return json.dumps(out, indent=4)

    def export(self, directory):
        """
        """
        os.chdir(directory)
        os.system("mkdir -p post-processing")
        os.chdir("post-processing")
        self.plot_info()
        self.markdown_report()
        if self.params["RUN_TYPE"].upper() == "GEO_OPT":
            self.get_optimized_structure_geo_opt(xyztraj="../ab-initio-pos-1.xyz", geo_opt_inp="../geo-opt.inp", directory="./")
            with open("geo-opt.json", 'w') as fout:
                fout.write(self.to_json_string())
        elif self.params["RUN_TYPE"].upper() == "CELL_OPT":
            self.get_optimized_structure_cell_opt(xyztraj="../ab-initio-pos-1.xyz", cell_opt_out="../cell-opt.out", directory="./")
            with open("cell-opt.json", 'w') as fout:
                fout.write(self.to_json_string())            
        os.chdir("../../")


class opt_post:
    """
    """
    def __init__(self, output, run_type):
        self.file = output
        self.run_type = run_type
        self.program_ended = None # whether calculation is finished
        self.cell = None # optimized cell
        self.atoms = None # optimized atoms
        self.opt_params = {}
        self.run_info = {}

        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()
        self.get_info()

    def get_info(self, output="geo-opt.out"):
        """
        get the general information of opt run from opt run output file
        which is now stored in self.lines
        """
        # check whether calculation has finished
        # cp2k will output the calculation directory on the last line
        # and different length of directory string will affect the content
        # of the last line, so we can't directly judge whether the calculation
        # has stopped only from the last line. and the solution is to check
        # the last several lines
        self.program_ended = False
        for line in self.lines[-9:]:
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "****" and line.split()[1] == "****" and line.split()[2] == "******" and line.split()[3] == "**" and line.split()[4] == "PROGRAM":
                self.program_ended = True
                break
        #
        self.get_final_structure()

        self.get_opt_params_and_run_info()

    def get_final_structure(self):
        pass

    def get_opt_params_and_run_info(self):
        """
        """
        self.run_info["scf-steps-converged"] = []
        self.run_info["total-energies"] = []
        self.run_info["fermi-energies"] = []

        for line in self.lines:
            # if it is an empty line continue to next line
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "****" and line.split()[1] == "****" and line.split()[5] == "STARTED":
                self.run_info["start-time"] = line.split("\n")[0]
            if line.split()[0] == "****" and line.split()[1] == "****" and line.split()[5] == "ENDED":
                self.run_info["stop-time"] = line.split("\n")[0]
            if line.split()[0] == "eps_scf:":
                self.opt_params["EPS_SCF"] = float(line.split()[1])
            if line.split()[0] == "Mixing" and line.split()[1] == 'method:':
                self.opt_params["MIXING"] = line.split()[2]
            if line.split()[0] == "added" and line.split()[1] == "MOs":
                self.opt_params["ADDED_MOS"] = int(line.split()[2])
            # there might be scf cycles that are not converged
            # so the total scf cycles might be more than the len(self.run_info["scf-steps-converged"])
            if line.split()[0] == "***" and line.split()[1] == "SCF" and line.split()[3] == "converged":
                self.run_info["scf-steps-converged"].append(int(line.split()[5]))
            # in every scf cycle, whether it converged or not, the ENERGY| Total FORCE_EVAL ( QS  ) energy (a.u.):
            # will always be print out
            # so len(self.run_info["total-energies"]) equal to the actually times for scf calculation
            if line.split()[0] == "ENERGY|" and line.split()[4] == "QS":
                self.run_info["total-energies"].append(float(line.split()[8]))

        #self.run_info["scf-cycles"] = len(self.run_info["scf-steps"])
        # the total scf cycles might be more than the len(self.run_info["scf-steps"])
        # so we should use self.run_info["total-energies"] to get the total number of scf cycles
        self.run_info["scf-cycles-total"] = len(self.run_info["total-energies"])

    def view_trajectory(self, trajfile="ab-initio-pos-1.xyz"):
        #os.system("xcrysden --xyz %s" % trajfile)
        subprocess.call(["xcrysden", "--xyz", trajfile])

    def plot_run_info(self):
        """
        """
        plt.plot(self.run_info["scf-steps-converged"])
        plt.title("Iterations per SCF(converged)")
        plt.xlabel("Scf cycle(converged)")
        plt.ylabel("Scf iterations")
        plt.tight_layout()
        plt.savefig("iterations-per-scf-converged.png")
        plt.close()

        plt.plot(self.run_info["total-energies"])
        plt.title("Total energies per SCF")
        plt.xlabel("Scf cycle")
        plt.ylabel("Total Energies (a.u.)")
        plt.tight_layout()
        plt.savefig("total-energies-per-scf.png")
        plt.close()

        plt.plot(self.run_info["fermi-energies"])
        plt.title("Fermi energies per SCF")
        plt.xlabel("Scf cycle")
        plt.ylabel("Fermi energies ()")
        plt.tight_layout()
        plt.savefig("fermi-energies-per-scf.png")
        plt.close()

    def markdown_report(self, md="OptimizationReport.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# 几何优化实验统计\n")
            fout.write("几何优化类型: %s\n" % self.run_type)
            fout.write("几何优化是否结束: %s\n" % self.program_ended)
            fout.write("## 优化参数\n")
            for item in self.opt_params:
                fout.write("- %s: %s\n" % (item, str(self.opt_params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out
            # note that the accuracy for the seconds is not fully guranteed 
            # e.g. 2019-11-26 12:09:36.487 is read as 2019-11-26 12:09:36
            start = datetime.datetime.strptime(self.run_info["start-time"].split()[7]+"-"+self.run_info["start-time"].split()[8].split(".")[0], "%Y-%m-%d-%H:%M:%S")
            if self.program_ended == True:
                stop = datetime.datetime.strptime(self.run_info["stop-time"].split()[7]+"-"+self.run_info["stop-time"].split()[8].split(".")[0], "%Y-%m-%d-%H:%M:%S")
                delta_t = stop -start
            fout.write("- Time consuming:\n")
            if self.program_ended == True:
                fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
            else:
                fout.write("  - job is not finished yet, but it starts at %s\n" % start)
            # end the time information
            for item in self.run_info:
                fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

            fout.write("## 运行信息图示\n")
            fout.write("Iterations per SCF(converged)\n")
            fout.write("![Iterations per SCF](iterations-per-scf-converged.png)\n")
            
            fout.write("Total energies per SCF\n")
            fout.write("![Total energies per SCF](total-energies-per-scf.png)\n")

            fout.write("Fermi energies per SCF\n")
            fout.write("![Fermi energies per SCF](fermi-energies-per-scf.png)\n")


    def export(self):
        self.plot_run_info()
        self.markdown_report("OptimizationReport.md")
