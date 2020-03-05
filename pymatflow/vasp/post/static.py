import datetime
import matplotlib.pyplot as plt



class scf_out:
    """
    """
    def __init__(self):
        """
        output is the output file of scf run
        """
        self.run_params = {}
        self.run_info = {}
        self.job_done = None


    def get_info(self, outcar):
        """
        get the general information of scf run from scf run output file
        which is now stored in self.lines
        """
        self.outcar = outcar
        with open(self.outcar, 'r') as fin:
            self.lines = fin.readlines()
        # check whether calculation is finished
        if len(self.lines[-1].split()) == 4 and self.lines[-1].split()[0] == "Voluntary" and self.lines[-1].split()[1] == "context":
            self.job_done = True
        else:
            self.job_done = False
        self.get_run_params_and_run_info()
    #
    def get_run_params_and_run_info(self):
        """
        """
        self.run_info["scf_energies"] = []

        for line in self.lines:
            # if it is an empty line continue to next line
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "executed" and line.split()[1] == "on" and line.split()[3] == "date":
                self.run_info["start_time"] = line.split("\n")[0]
            #if line.split()[0] == "This" and line.split()[1] == "run" and line.split()[3] == "terminated":
            #    self.run_info["stop-time"] = line.split("\n")[0]
            if line.split()[0] == "Total" and line.split()[1] == "CPU" and line.split()[2] == "time":
                self.run_info["total_cpu_time"] = float(line.split()[5]) # in unit of second
            if line.split()[0] == "Elapsed" and line.split()[1] == "time":
                self.run_info["elapsed_time"] = float(line.split()[3])
            if line.split()[0] == "energy" and line.split()[1] == "without" and line.split()[2] == "entropy":
                self.run_info["scf_energies"].append(float(line.split()[4]))
            if line.split()[0] ==  "E_fermi" and line.split()[1] == ":":
                self.run_info["fermi_energy"] = float(line.split()[2])
            if line.split()[0] == "FORCES:" and line.split()[1] == "max":
                self.run_info["forces_rms"] = float(line.split()[5])
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




class scf_post:
    """
    """
    def __init__(self):
        """
        output is the output file of scf run
        """
        self.electronic_params = {}
        self.ionic_params = {}
        self.scf_params = {}
        self.run_info = {}
        self.job_done = None

    def get_outcar(self, outcar):
        self.outcar = outcar
        with open(self.outcar, 'r') as fin:
            self.lines = fin.readlines()
        self.get_info()

    def get_info(self):
        """
        get the general information of scf run from scf run output file
        which is now stored in self.lines
        """
        # check whether calculation is finished
        if len(self.lines[-1].split()) == 4 and self.lines[-1].split()[0] == "Voluntary" and self.lines[-1].split()[1] == "context":
            self.job_done = True
        else:
            self.job_done = False
        self.get_scf_params_and_run_info()
    #
    def get_scf_params_and_run_info(self):
        """
        """
        self.run_info["scf-energies"] = []

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
            if line.split()[0] == "energy" and line.split()[1] == "without" and line.split()[2] == "entropy":
                self.run_info["scf-energies"].append(float(line.split()[4]))
            if line.split()[0] ==  "E-fermi" and line.split()[1] == ":":
                self.run_info["fermi-energy"] = float(line.split()[2])
            if line.split()[0] == "FORCES:" and line.split()[1] == "max":
                self.run_info["forces-rms"] = float(line.split()[5])
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

    def plot_run_info(self):
        """
        """
        plt.plot(self.run_info["scf-energies"])
        plt.title("Energy per scf step")
        plt.xlabel("Scf step")
        plt.ylabel("Total energy")
        plt.tight_layout()
        plt.savefig("energy-per-scf-step.png")
        plt.close()



    def markdown_report(self, md="SCFReport.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# SCF实验统计\n")
            fout.write("## SCF参数\n")
            for item in self.scf_params:
                fout.write("- %s: %s\n" % (item, str(self.scf_params[item])))
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
            fout.write("Total energy per scf step\n")
            fout.write("![energy per scf step](energy-per-scf-step.png)\n")

    def export(self):
        self.plot_run_info()
        self.markdown_report("SCFReport.md")
