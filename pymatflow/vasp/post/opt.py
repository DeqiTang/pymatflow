import os
import sys
import datetime
import subprocess
import matplotlib.pyplot as plt

from pymatflow.vasp.base.xyz import Atom

class opt_out:
    """
    Note:
    """
    def __init__(self):
        """
        output:
        """
        self.outcar = None
        self.poscar = None
        self.run_params = {}
        self.run_info = {}
        self.job_done = None # whether calculation has finished
        self.relaxed = None # whether structure is relaxed successfully

        self.cell = None # optimized cell
        self.trajectory = None



    def get_info(self, outcar, poscar):
        """
        get the general information of opt run from opt run output file
        which is now stored in self.lines
        """
        self.outcar = outcar
        self.poscar = poscar
        with open(self.outcar, 'r') as fout:
            self.lines = fout.readlines()
        with open(self.poscar, 'r') as fin:
            self.lines_poscar = fin.readlines()
        # check whether calculation is finished
        if len(self.lines[-1].split()) == 4 and self.lines[-1].split()[0] == "Voluntary" and self.lines[-1].split()[1] == "context":
            self.job_done = True
        else:
            self.job_done = False
        # check whether successfully relaxed
        self.relaxed = False
        for line in self.lines:
            if line == " reached required accuracy - stopping structural energy minimisation\n":
                self.relaxed = True
                break

        self.get_trajectory()
        self.get_opt_params_and_run_info()


    def get_trajectory(self):
        """
        in OUTCAR we can get the atom positions but the element name is not available
        so we have to found it in POSCAR through self.lines_poscar

        Note: the cell of each trajectory is also obtained
        """

        self.trajectory = []
        self.cells = []
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "POSITION" and self.lines[i].split()[1] == "TOTAL-FORCE":
                atm = []
                j = i + 2
                while len(self.lines[j].split()) == 6:
                    atm.append(Atom("x", float(self.lines[j].split()[0]), float(self.lines[j].split()[1]), float(self.lines[j].split()[2])))
                    j = j + 1
                self.trajectory.append(atm)
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "direct" and self.lines[i].split()[1] == "lattice" and len(self.lines[i-4].split()) > 0 and self.lines[i-4].split()[0] == "VOLUME":
                # note in the begining vasp will output the direct lattice vector of the input structure so that must be excluded, so we check the previous 4 line by i-4 above
                cell = []
                for j in range(3):
                    cell.append([float(self.lines[i+1].split()[0]), float(self.lines[i+1].split()[1]), float(self.lines[i+1].split()[2])])
                    cell.append([float(self.lines[i+2].split()[0]), float(self.lines[i+2].split()[1]), float(self.lines[i+2].split()[2])])
                    cell.append([float(self.lines[i+3].split()[0]), float(self.lines[i+3].split()[1]), float(self.lines[i+3].split()[2])])
                self.cells.append(cell)
        # len(self.cells) must equals to len(self.trajectory)
        if not len(self.cells) == len(self.trajectory):
            print(len(self.cells))
            print(len(self.trajectory))
            print("============================================================\n")
            print("                    WARNING!!!!!!!\n")
            print("------------------------------------------------------------\n")
            print("post.opt.opt_out.get_trajectory():\n")
            print("length of self.cells and self.trajectory should be equal!!!!\n")
            print('but they are not now! please check it manually\n')
            sys.exit(1)
        # set labels for each atom in each image
        elements = self.lines_poscar[5].split() # like ["S", "Mo"]
        natom_each_element = self.lines_poscar[6].split() # like ["2", "1"]
        label_list = []
        for i in range(len(elements)):
            for j in range(int(natom_each_element[i])):
                label_list.append(elements[i])
        for i in range(len(self.trajectory)):
            for j in range(len(label_list)):
                self.trajectory[i][j].name = label_list[j]
        #


    def get_opt_params_and_run_info(self):
        """
        run_info["iterations"]: scf iterations per scf step
        run_info["total_energies"]: total energies of every scf step
        run_info["fermi_energies"]: fermi energies of every scf step
        run_info["total_forces_rms"]: total RMS forces of every scf step
        """
        self.run_info["iterations"] = []
        self.run_info["total_energies"] = []
        self.run_info["fermi_energies"] = []
        self.run_info["total_forces_rms"] = []

        for line in self.lines:
            # if it is an empty line continue to next line
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "executed" and line.split()[1] == "on" and line.split()[3] == "date":
                self.run_info["start_time"] = line.split("\n")[0]
            #if line.split()[0] == "This" and line.split()[1] == "run" and line.split()[3] == "terminated":
            #    self.run_info["stop_time"] = line.split("\n")[0]
            if line.split()[0] == "Total" and line.split()[1] == "CPU" and line.split()[2] == "time":
                self.run_info["total_cpu_time"] = float(line.split()[5]) # in unit of second
            if line.split()[0] == "Elapsed" and line.split()[1] == "time":
                self.run_info["elapsed_time"] = float(line.split()[3])
            if len(line.split()) == 4 and line.split()[1] == "Iteration":
                self.run_info["iterations"].append(line)
            if line.split()[0] == "energy" and line.split()[1] == "without" and line.split()[2] == "entropy=":
                self.run_info["total_energies"].append(float(line.split()[3]))
            if line.split()[0] ==  "E-fermi" and line.split()[1] == ":":
                self.run_info["fermi_energies"].append(float(line.split()[2]))
            if line.split()[0] == "FORCES:" and line.split()[1] == "max":
                self.run_info["total_forces_rms"].append(float(line.split()[5]))
            if line.split()[0] == "ENCUT" and line.split()[1] == "=":
                self.run_params["ENCUT"] = float(line.split()[2])
            if line.split()[0] == "EDIFF" and line.split()[1] == "=":
                self.run_params["EDIFF"] = float(line.split()[2])
            if line.split()[0] == "LREAL" and line.split()[1] == "=":
                self.run_params["LREAL"] = line.split()[2]
            if line.split()[0] == "EDIFFG" and line.split()[1] == "=":
                self.run_params["EDIFFG"] = float(line.split()[2])
            if line.split()[0] == "NSW" and line.split()[1] == "=":
                self.run_params["NSW"] = int(line.split()[2])
            if line.split()[0] == "IBRION" and line.split()[1] == "=":
                self.run_params["IBRION"] = int(line.split()[2])
            if line.split()[0] == "NFREE" and line.split()[1] == "=":
                self.run_params["NFREE"] = int(line.split()[2])
            if line.split()[0] == "ISIF" and line.split()[1] == "=":
                self.run_params["ISIF"] = int(line.split()[2])
            if line.split()[0] == "POTIM" and line.split()[1] == "=":
                self.run_params["POTIM"] = float(line.split()[2])
            if line.split()[0] == "TEIN" and line.split()[1] == "=":
                self.run_params["TEIN"] = float(line.split()[2])
            if line.split()[0] == "TEBEG" and line.split()[1] == "=":
                self.run_params["TEBEG"] = float(line.split()[2].split(";")[0])
            if line.split()[0] == "SMASS" and line.split()[1] == "=":
                self.run_params["SMASS"] = float(line.split()[2])
            if line.split()[0] == "PSTRESS=":
                self.run_params["PSTRESS"] = float(line.split()[1])

        #self.run_info["scf-cycles"] = len(self.run_info["iterations"])
        #if self.run_type == "relax":
        #    self.run_info["ion-steps"] = len(self.run_info["iterations"]) - 1
        #elif self.run_type == "vc-relax":
        #    self.run_info["ion-steps"] = len(self.run_info["iterations"]) - 2

    def print_final_structure(self, xyz="optimized.xyz"):
        if self.relaxed == False:
            with open("final-structure(not-relaxed).xyz", 'w') as fout:
                fout.write("%d\n" % len(self.trajectory[0]))
                fout.write("Warning: structure failed to be relaxed or vc-relaxed, this is the final structure(unrelaxed)\n")
                for atom in self.trajectory[-1]:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
            return
        # printout relaxed structure
        with open(xyz, 'w') as fout:
            fout.write("%d\n" % len(self.trajectory[-1]))
            fout.write("cell: %f %f %f | %f %f %f | %f %f %f\n" % (
                self.cells[-1][0][0], self.cells[-1][0][1], self.cells[-1][0][2],
                self.cells[-1][1][0], self.cells[-1][1][1], self.cells[-1][1][2],
                self.cells[-1][2][0], self.cells[-1][2][1], self.cells[-1][2][2]
                ))
            for atom in self.trajectory[-1]:
                fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))

    def print_trajectory(self, xyz="trajectory.xyz"):
        with open(xyz, 'w') as fout:
            for i in range(len(self.trajectory)):
                fout.write("%d\n" % len(self.trajectory[i]))
                #fout.write("i = %d\n" % i)
                fout.write("cell(i=%d): %f %f %f | %f %f %f | %f %f %f\n" % (
                    i,
                    self.cells[i][0][0], self.cells[i][0][1], self.cells[i][0][2],
                    self.cells[i][1][0], self.cells[i][1][1], self.cells[i][1][2],
                    self.cells[i][2][0], self.cells[i][2][1], self.cells[i][2][2]
                    ))
                for atom in self.trajectory[i]:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))

    def export(self, directory="tmp-vasp-optimization"):
        os.chdir(directory)
        os.system("mkdir -p post-processing")
        os.chdir("post-processing")
        # now we are in post-processing, generate the output and return
        self.print_trajectory()
        self.print_final_structure()

        #plt.plot(self.run_info["self-energies"])
        #plt.title("Energy per scf step")
        #plt.xlabel("Scf step")
        #plt.ylabel("Total energy")
        #plt.tight_layout()
        #plt.savefig("energy-per-scf-step.png")
        #plt.close()

        with open("self-info.md", 'w', encoding='utf-8') as fout:
            fout.write("# 几何优化实验统计\n")
            fout.write("几何优化类型: ISIF = %d\n" % self.run_params["ISIF"])
            fout.write("几何优化任务是否结束:%s\n" % str(self.job_done))
            if self.job_done == True:
                fout.write("是否成功优化: %s\n" % str(self.relaxed))
            else:
                fout.write("是否成功优化: %s\n" % ("运行未结束, 结果未知"))
            fout.write("## 离子步参数\n")
            for item in self.run_params:
                fout.write("- %s: %s\n" % (item, str(self.run_params[item])))
            fout.write("## 电子步参数\n")
            for item in self.run_params:
                fout.write("- %s: %s\n" % (item, str(self.run_params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out
            # Importante: the length of the time string might be different, depending
            # on the value of hours and minutes and seconds. if they are two digits
            # number, they will be divided like: '11: 6: 2', only when they all are
            # two digtis number, they will not be divided '11:16:12'
            # so we have to preprocess it to build the right time string to pass into
            # datetime.datetime.strptime()
            start_str = self.run_info["start_time"].split()[4]+"-"+self.run_info["start_time"].split()[5]
            if self.job_done == True:
                #stop_str = self.run_info["stop-time"].split()[8]+"-"+self.run_info["stop-time"].split()[5]+self.run_info["stop-time"].split()[6]+self.run_info["stop-time"].split()[7]
                pass

            start = datetime.datetime.strptime(start_str, "%Y.%m.%d-%H:%M:%S")
            #if self.job_done == True:
            #    stop = datetime.datetime.strptime(stop_str, "%d%b%Y-%H:%M:%S")
            #    delta_t = stop -start
            fout.write("- Time consuming:\n")
            fout.write("  - job starts at %s\n" % start)
            fout.write("  - Elapsed time: %.3f(sec) = %.3f(min) = %.3f(hour)\n" % (self.run_info["elapsed_time"], self.run_info["elapsed_time"]/60, self.run_info["elapsed_time"]/3600))
            #if self.job_done == True:
            #    fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
            #else:
            #    fout.write("  - job is not finished yet, but it starts at %s\n" % start)
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
            fout.write("![Total forces per SCF](total-forces-rms-per-scf.png)\n")

        os.chdir("../")
        os.chdir("../")


class opt_post:
    """
    Note:
    """
    def __init__(self, output="OUTCAR"):
        """
        output:
        """
        self.file = output
        self.electronic_params = {}
        self.ionic_params = {}
        self.run_info = {}
        self.job_done = None # whether calculation has finished
        self.relaxed = None # whether structure is relaxed successfully

        self.cell = None # optimized cell
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
        if len(self.lines[-1].split()) == 4 and self.lines[-1].split()[0] == "Voluntary" and self.lines[-1].split()[1] == "context":
            self.job_done = True
        else:
            self.job_done = False
        # check whether successfully relaxed
        self.relaxed = False
        for line in self.lines:
            if line == " reached required accuracy - stopping structural energy minimisation\n":
                self.relaxed = True
                break

        self.get_trajectory()
        self.get_opt_params_and_run_info()


    def get_trajectory(self):
        self.trajectory = []
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "POSITION" and self.lines[i].split()[1] == "TOTAL-FORCE":
                atm = []
                j = i + 2
                while len(self.lines[j].split()) == 6:
                    atm.append(Atom("x", float(self.lines[j].split()[0]), float(self.lines[j].split()[1]), float(self.lines[j].split()[2])))
                    j = j + 1
                self.trajectory.append(atm)
        #
        """
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
        """

    def get_opt_params_and_run_info(self):
        """
        run_info["iterations"]: scf iterations per scf step
        run_info["total-energies"]: total energies of every scf step
        run_info["fermi-energies"]: fermi energies of every scf step
        run_info["total-forces-rms"]: total RMS forces of every scf step
        """
        self.run_info["iterations"] = []
        self.run_info["total-energies"] = []
        self.run_info["fermi-energies"] = []
        self.run_info["total-forces-rms"] = []

        for line in self.lines:
            # if it is an empty line continue to next line
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "executed" and line.split()[1] == "on" and line.split()[3] == "date":
                self.run_info["start-time"] = line.split("\n")[0]
            #if line.split()[0] == "This" and line.split()[1] == "run" and line.split()[3] == "terminated":
            #    self.run_info["stop-time"] = line.split("\n")[0]
            if line.split()[0] == "Total" and line.split()[1] == "CPU" and line.split()[2] == "time":
                self.run_info["total-cpu-time"] = float(line.split()[5]) # in unit of second
            if line.split()[0] == "Elapsed" and line.split()[1] == "time":
                self.run_info["elapsed-time"] = float(line.split()[3])
            if len(line.split()) == 4 and line.split()[1] == "Iteration":
                self.run_info["iterations"].append(line)
            if line.split()[0] == "energy" and line.split()[1] == "without" and line.split()[2] == "entropy=":
                self.run_info["total-energies"].append(float(line.split()[3]))
            if line.split()[0] ==  "E-fermi" and line.split()[1] == ":":
                self.run_info["fermi-energies"].append(float(line.split()[2]))
            if line.split()[0] == "FORCES:" and line.split()[1] == "max":
                self.run_info["total-forces-rms"].append(float(line.split()[5]))
            if line.split()[0] == "ENCUT" and line.split()[1] == "=":
                self.electronic_params["ENCUT"] = float(line.split()[2])
            if line.split()[0] == "EDIFF" and line.split()[1] == "=":
                self.electronic_params["EDIFF"] = float(line.split()[2])
            if line.split()[0] == "LREAL" and line.split()[1] == "=":
                self.electronic_params["LREAL"] = line.split()[2]
            if line.split()[0] == "EDIFFG" and line.split()[1] == "=":
                self.ionic_params["EDIFFG"] = float(line.split()[2])
            if line.split()[0] == "NSW" and line.split()[1] == "=":
                self.ionic_params["NSW"] = int(line.split()[2])
            if line.split()[0] == "IBRION" and line.split()[1] == "=":
                self.ionic_params["IBRION"] = int(line.split()[2])
            if line.split()[0] == "NFREE" and line.split()[1] == "=":
                self.ionic_params["NFREE"] = int(line.split()[2])
            if line.split()[0] == "ISIF" and line.split()[1] == "=":
                self.ionic_params["ISIF"] = int(line.split()[2])
            if line.split()[0] == "POTIM" and line.split()[1] == "=":
                self.ionic_params["POTIM"] = float(line.split()[2])
            if line.split()[0] == "TEIN" and line.split()[1] == "=":
                self.ionic_params["TEIN"] = float(line.split()[2])
            if line.split()[0] == "TEBEG" and line.split()[1] == "=":
                self.ionic_params["TEBEG"] = float(line.split()[2].split(";")[0])
            if line.split()[0] == "SMASS" and line.split()[1] == "=":
                self.ionic_params["SMASS"] = float(line.split()[2])
            if line.split()[0] == "PSTRESS=":
                self.ionic_params["PSTRESS"] = float(line.split()[1])

        #self.run_info["scf-cycles"] = len(self.run_info["iterations"])
        #if self.run_type == "relax":
        #    self.run_info["ion-steps"] = len(self.run_info["iterations"]) - 1
        #elif self.run_type == "vc-relax":
        #    self.run_info["ion-steps"] = len(self.run_info["iterations"]) - 2

    def print_final_structure(self, xyz="optimized.xyz"):
        if self.relaxed == False:
            with open("final-structure(not-relaxed).xyz", 'w') as fout:
                fout.write("%d\n" % len(self.trajectory[0]))
                fout.write("Warning: structure failed to be relaxed or vc-relaxed, this is the final structure(unrelaxed)\n")
                for atom in self.trajectory[-1]:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
            return
        # printout relaxed structure
        with open(xyz, 'w') as fout:
            fout.write("%d\n" % len(self.trajectory[-1]))
            fout.write("cell:\n")
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
        #plt.plot(self.run_info["iterations"])
        #plt.title("Iterations per SCF")
        #plt.xlabel("Scf cycles")
        #plt.ylabel("iterations")
        #plt.tight_layout()
        #plt.savefig("iterations-per-scf.png")
        #plt.close()

        plt.plot(self.run_info["total-energies"])
        plt.title("Total energies per SCF")
        plt.xlabel("Scf cycles")
        plt.ylabel("Total Energies (eV)")
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

        plt.plot(self.run_info["total-forces-rms"])
        plt.title("Total forces(RMS) per SCF")
        plt.xlabel("Scf cycles")
        plt.ylabel("Total forces (eV/Angst)")
        plt.tight_layout()
        plt.savefig("total-forces-rms-per-scf.png")
        plt.close()

    def markdown_report(self, md="OptimizationReport.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# 几何优化实验统计\n")
            fout.write("几何优化类型: ISIF = %d\n" % self.ionic_params["ISIF"])
            fout.write("几何优化任务是否结束:%s\n" % str(self.job_done))
            if self.job_done == True:
                fout.write("是否成功优化: %s\n" % str(self.relaxed))
            else:
                fout.write("是否成功优化: %s\n" % ("运行未结束, 结果未知"))
            fout.write("## 离子步参数\n")
            for item in self.ionic_params:
                fout.write("- %s: %s\n" % (item, str(self.ionic_params[item])))
            fout.write("## 电子步参数\n")
            for item in self.electronic_params:
                fout.write("- %s: %s\n" % (item, str(self.electronic_params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out
            # Importante: the length of the time string might be different, depending
            # on the value of hours and minutes and seconds. if they are two digits
            # number, they will be divided like: '11: 6: 2', only when they all are
            # two digtis number, they will not be divided '11:16:12'
            # so we have to preprocess it to build the right time string to pass into
            # datetime.datetime.strptime()
            start_str = self.run_info["start-time"].split()[4]+"-"+self.run_info["start-time"].split()[5]
            if self.job_done == True:
                #stop_str = self.run_info["stop-time"].split()[8]+"-"+self.run_info["stop-time"].split()[5]+self.run_info["stop-time"].split()[6]+self.run_info["stop-time"].split()[7]
                pass

            start = datetime.datetime.strptime(start_str, "%Y.%m.%d-%H:%M:%S")
            #if self.job_done == True:
            #    stop = datetime.datetime.strptime(stop_str, "%d%b%Y-%H:%M:%S")
            #    delta_t = stop -start
            fout.write("- Time consuming:\n")
            fout.write("  - job starts at %s\n" % start)
            fout.write("  - Elapsed time: %.3f(sec) = %.3f(min) = %.3f(hour)\n" % (self.run_info["elapsed-time"], self.run_info["elapsed-time"]/60, self.run_info["elapsed-time"]/3600))
            #if self.job_done == True:
            #    fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
            #else:
            #    fout.write("  - job is not finished yet, but it starts at %s\n" % start)
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
            fout.write("![Total forces per SCF](total-forces-rms-per-scf.png)\n")


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
