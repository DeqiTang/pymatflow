#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys

import numpy as np
import matplotlib.pyplot as plt



class pdos_out:
    """
    """
    def __init__(self):
        self.data = {} # contain the pdos data but not the tdos
        self.tdos = None
        self.energies = None

    def get_data(self, directory="tmp-qe-static", filpdos="projwfc"):
        """
        this function first try to get fermi energy from the nscfout file
        and if nscfout doesn't exist it will try to extract fermi energy
        from scfout. if both don't exist, it will stop the program and
        print out the warnnig, which guarantee that the fermi energy is
        always shifted to 0

        atomorb is a string in format like this atm#1(Li)_wfc#2(s).
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("pdos post:\n")
            print("  directory of previous scf calculattion not found!\n")
            sys.exit(1)

        os.chdir(directory)
        os.system("ls | grep '%s.pdos_' > projwfc-pdos-file.data" % filpdos)
        with open("projwfc-pdos-file.data", 'r') as fin:
            for line in fin:
                if line.split(".")[1] == "pdos_tot\n":
                    with open(line.split()[0], 'r') as f:
                        f.readline()
                        self.tdos = np.loadtxt(f)
                    continue
                atmorb = line.split("_")[1]+"_"+line.split("_")[2].split()[0]
                # atomorb is a string in format like this atm#1(Li)_wfc#2(s).
                with open(line.split()[0], 'r') as f:
                    f.readline()
                    self.data[atmorb] = np.loadtxt(f)
        # check information on spin
        with open("%s.pdos_tot" % filpdos, 'r') as fin:
            first_line = fin.readline()
        if "pdosup(E)" in first_line.split() and "pdosdw(E)" in first_line.split():
            if "dosup(E)" in first_line.split() and "dosdw(E)" in first_line.split():
                self.magnetic_status = "collinear-spin-polarized"
            else:
                self.magnetic_status = "non-collinear-non-spin-orbit"
        else:
            self.magnetic_status = "collinear-spin-unpolarized"
        os.chdir("../")

        self.energies = self.tdos[:, 0]

        # get fermi energy from nscf output
        scfout = "static-scf.out"
        nscfout = "static-nscf.out"
        if os.path.exists(os.path.join(directory, nscfout)):
            with open(os.path.join(directory, nscfout), 'r') as fin:
                for line in fin:
                    if len(line.split()) == 0:
                        continue
                    if line.split()[0] == "the" and line.split()[1] == "Fermi":
                        efermi = float(line.split()[4])
        elif os.path.exists(os.path.join(directory, scfout)):
            with open(os.path.join(directory, scfout), 'r') as fin:
                for line in fin:
                    if len(line.split()) == 0:
                        continue
                    if line.split()[0] == "the" and line.split()[1] == "Fermi":
                        efermi = float(line.split()[4])
        else:
            print("===========================================================\n")
            print("                Warning !!!\n")
            print("===========================================================\n")
            print("PDOS postprocessing:\n")
            print("must provide nscfout or at least scfout to get Fermi energy\n")
            sys.exit(1)
        # shift fermie energy to 0
        for i in range(len(self.energies)):
            self.energies[i] = self.energies[i] - efermi
        print("===============================================\n")
        print("qe.post.pdos:\n")
        print("we automatically shift the fermi energy\n")
        print("from %f to 0\n" % efermi)
        print("efermi is read from static-nscf.out\n")
        print("or statis-scf.out, if static-nscf.out is not available\n")
        #
        # tranfser self.data to new data structure for better usage:
        self._transfer_data_to_new_struct()

    def _transfer_data_to_new_struct(self):
        """
        self.magnetic_status:
            "collinear-spin-unpolarized" -> self.data_0
            "collinear-spin-polarized" -> self.data_1
            "non-collinear-non-spin-orbit" -> self.data_2
            "non-collinear-spin-orbit" -> self.data_3
        """
        if self.magnetic_status == "collinear-spin-unpolarized":
            """
            # pdos are in format like this:
            # E LDOS(E) PDOS_1(E) ... PDOS_2l+1(E)
            # self.data_0:
            {
                atmorb: {
                    "ldos": [],
                    "pdos_l_1": [],
                    "pdos_l_2": [],
                    ...
                    "pdos_l_2l+1": []
                }
            }
            l could be: s, p, d, f, in format like this: 1(s), 2(s), 2(p)
            self.data_0_tdos:
            {
                dos: [],
                pdos: [],
            }
            """
            self.data_0 = {}
            for atmorb in self.data:
                self.data_0[atmorb] = {}
                self.data_0[atmorb]["ldos"] = self.data[atmorb][:, 1]
                for i in range(self.data[atmorb].shape[1]-2):
                    self.data_0[atmorb]["pdos_%s_%d" % (self.get_orb_type(atmorb), i+1)] = self.data[atmorb][:, i+2]
            #
            self.data_0_tdos = {}
            self.data_0_tdos["dos"] = self.tdos[:, 1]
            self.data_0_tdos["pdos"] = self.tdos[:, 2]

        elif self.magnetic_status == "collinear-spin-polarized":
            """
            # pdos are in format like this:
            # E ldosup(E) ldosdw(E) pdos_1up(E) pdos_1dw(E) ... pdow_2l+1up(E) pdos_2l+1dw(E)
            # self.data_1:
            {
                atmorb: {
                    "ldos_up": [],
                    "ldos_dw": [],
                    "pdos_l_1_up": [],
                    "pdos_l_1_dw": [],
                    "pdos_l_2_up": [],
                    "pdos_l_2_dw": [],
                    ...
                    "pdos_l_2l+1_up": [],
                    "pdos_l_2l+1_dw": [],
                }
            }
            l could be: s, p, d, f, in format like this: 1(s), 2(s), 2(p)
            self.data_1_tdos:
            {
                "dos_up": [],
                "dos_dw": [],
                "pdos_up": [],
                "pdos_dw": []
            }
            """
            self.data_1 = {}
            for atmorb in self.data:
                self.data_1[atmorb] = {}
                self.data_1[atmorb]["ldos_up"] = self.data[atmorb][:, 1]
                self.data_1[atmorb]["ldos_dw"] = self.data[atmorb][:, 2]
                for i in range(int((self.data[atmorb].shape[1]-3)/2)):
                    self.data_1[atmorb]["pdos_%s_%d_up" % (self.get_orb_type(atmorb), i+1)] = self.data[atmorb][:, 3+2*i]
                    self.data_1[atmorb]["pdos_%s_%d_dw" % (self.get_orb_type(atmorb), i+1)] = self.data[atmorb][:, 3+2*i+1]
            #
            self.data_1_tdos = {}
            self.data_1_tdos["dos_up"] = self.tdos[:, 1]
            self.data_1_tdos["dos_dw"] = self.tdos[:, 2]
            self.data_1_tdos["pdos_up"] = self.tdos[:, 3]
            self.data_1_tdos["pdos_dw"] = self.tdos[:, 4]
        elif self.magnetic_status == "non-collinear-non-spin-orbit":
            pass
        elif self.magnetic_status == "non-collinear-spin-orbit":
            pass


    def plot_elem_orb_l_proj(self, plotrange=[0.0, 1.0], filename="pdos-projected-to-element-and-orbital-l.png", fontsize=10):
        """
        plotrange:
            a list of two values(between 0 and 1) defining the percentage
            of data to plot.
            plotrange[0]: left boundary of the data to plot
            plotrange[1]: right boundary of the data to plot
            default is plotrange[0] = 0, plotrange[1], in which case
            all the data will be plot.
        """
        if self.magnetic_status == "collinear-spin-unpolarized":
            data = {}
            for atmorb in self.data_0:
                key = self.get_elem_type(atmorb)+"-"+self.get_orb_type(atmorb)
                if key in data:
                    data[key] = data[key] + self.data_0[atmorb]["ldos"]
                else:
                    data[key] = self.data_0[atmorb]["ldos"]
        elif self.magnetic_status == "collinear-spin-polarized":
            data = {}
            for atmorb in self.data_1:
                key = self.get_elem_type(atmorb)+"-"+self.get_orb_type(atmorb)
                key_up = self.get_elem_type(atmorb)+"-"+self.get_orb_type(atmorb)+"-up"
                key_dw = self.get_elem_type(atmorb)+"-"+self.get_orb_type(atmorb)+"-down"
                if key_up in data and key_dw in data:
                    data[key_up] = data[key_up] + self.data_1[atmorb]["ldos_up"]
                    data[key_dw] = data[key_dw] + (-self.data_1[atmorb]["ldos_dw"])
                else:
                    data[key_up] = self.data_1[atmorb]["ldos_up"]
                    data[key_dw] = (-self.data_1[atmorb]["ldos_dw"])

        # plot the pdos in the specified percentage range
        begin = int(len(self.energies)*plotrange[0])
        end = int(len(self.energies)*plotrange[1])
        for key in data:
            plt.plot(self.energies[begin:end], data[key][begin:end], label=key)

        # plot the total dos in the specified percentage range
        if self.magnetic_status == "collinear-spin-unpolarized":
            plt.plot(self.energies[begin:end], self.data_0_tdos["dos"][begin:end], label="Total-DOS")
        elif self.magnetic_status == "collinear-spin-polarized":
            plt.plot(self.energies[begin:end], self.data_1_tdos["dos_up"], label="Total-DOS-Up")
            plt.plot(self.energies[begin:end], -self.data_1_tdos["dos_dw"], label="Total-DOS-Down")

        # set formats
        font = {'size': fontsize}
        plt.tick_params(labelsize=fontsize)

        plt.grid(which="major", axis="x", linewidth=0.75, linestyle="-", color="0.75")
        plt.grid(which="major", axis="y", linewidth=0.75, linestyle="-", color="0.75")
        plt.title("Projected Density of States")
        plt.xlabel(r"$\mathit{E}-\mathit{E}_\mathrm{f} \mathrm{(eV)}$", font)
        plt.ylabel("States", font)
        plt.legend(prop=font)
        plt.tight_layout()
        plt.savefig("%s" % filename)
        plt.close()


    def plot_atom_orb_l_proj(self, atomtoproj=[], plotrange=[0.0, 1.0], filename="pdos-projected-to-atom-and-orbital-l.png", fontsize=10):
        """
        plotrange:
            a list of two values(between 0 and 1) defining the percentage
            of data to plot.
            plotrange[0]: left boundary of the data to plot
            plotrange[1]: right boundary of the data to plot
            default is plotrange[0] = 0, plotrange[1], in which case
            all the data will be plot.
        atomtoproj:
            the list of atoms to do the projection. atom number starts with 1
        """
        # plot the data in the specified percentage range
        begin = int(len(self.energies)*plotrange[0])
        end = int(len(self.energies)*plotrange[1])

        # atom projected dos
        if self.magnetic_status == "collinear-spin-unpolarized":
            for atmorb in self.data_0:
                if self.get_atom_num(atmorb) in atomtoproj:
                    plt.plit(self.energies[begin:end], self.data_0[atmorb]["ldos"][begin:end], label="Atom(%d):%s-%s" % (self.get_atom_num(atmorb), self.get_elem_type(atmorb), self.get_orb_type(atmorb)))
        elif self.magnetic_status == "collinear-spin-polarized":
            for atmorb in self.data_1:
                if self.get_atom_num(atmorb) in atomtoproj:
                    plt.plot(self.energies[begin:end], self.data_1[atmorb]["ldos_up"][begin:end], label="Atom(%d):%s-%s-up" % (self.get_atom_num(atmorb), self.get_elem_type(atmorb), self.get_orb_type(atmorb)))
                    plt.plot(self.energies[begin:end], -self.data_1[atmorb]["ldos_dw"][begin:end], label="Atom(%d):%s-%s-down" % (self.get_atom_num(atmorb), self.get_elem_type(atmorb), self.get_orb_type(atmorb)))

        # set formats
        font = {'size': fontsize}
        plt.tick_params(labelsize=fontsize)

        plt.grid(which="major", axis="x", linewidth=0.75, linestyle="-", color="0.75")
        plt.grid(which="major", axis="y", linewidth=0.75, linestyle="-", color="0.75")
        plt.title("Projected(Atom) Density of States")
        plt.xlabel(r"$\mathit{E}-\mathit{E}_\mathrm{f} \mathrm{(eV)}$", font)
        plt.ylabel("States", font)
        if len(atomtoproj) != 0:
            plt.legend(prop=font)
        plt.tight_layout()
        plt.savefig("%s" % filename)
        plt.close()


    def plot_tdos(self, plotrange=[0, 1.0], filename="total-dos.png"):
        """
        plotrange:
            a list of two values(between 0 and 1) defining the percentage
            of data to plot.
            plotrange[0]: left boundary of the data to plot
            plotrange[1]: right boundary of the data to plot
            default is plotrange[0] = 0, plotrange[1], in which case
            all the data will be plot.
        """
        # plot the total dos in the specified percentage range
        begin = int(len(self.energies)*plotrange[0])
        end = int(len(self.energies)*plotrange[1])

        if self.magnetic_status == "collinear-spin-unpolarized":
            plt.plot(self.energies[begin:end], self.data_0_tdos["dos"][begin:end], label="Total-DOS")
        elif self.magnetic_status == "collinear-spin-polarized":
            plt.plot(self.energies[begin:end], self.data_1_tdos["dos_up"], label="Total-DOS-Up")
            plt.plot(self.energies[begin:end], -self.data_1_tdos["dos_dw"], label="Total-DOS-Down")

        plt.grid(which="major", axis="x", linewidth=0.75, linestyle="-", color="0.75")
        plt.grid(which="major", axis="y", linewidth=0.75, linestyle="-", color="0.75")
        plt.title("Total Density of States")
        plt.xlabel(r"$\mathit{E}-\mathit{E}_\mathrm{f} (eV)$")
        plt.ylabel("States")
        plt.legend()
        plt.tight_layout()
        plt.savefig("%s" % filename)
        plt.close()

    def get_elem_type(self, atmorb):
        """
        get element name from atmorb
        atmorb is the key in self.data
        it's like this:
            atm#1(Li)_wfc#2(s)
        return value of the above input
        will be:
            'Li'
        """
        return atmorb.split("(")[1].split(")")[0]

    def get_orb_type(self, atmorb):
        """
        get element name and orb from atmorb
        atmorb is the key in self.data
        it's like this:
            atm#1(Li)_wfc#2(s)
        return value of the above input
        will be:
            '2(s)'
        """
        return atmorb.split("#")[2]

    def get_atom_num(self, atmorb):
        """
        get atom name from atmorb
        atmorb is the key in self.data
        it's like this:
            atm#1(Li)_wfc#2(s)
        return value of the above input
        will be:
            1
        """
        return int(atmorb.split("(")[0].split("#")[1])

    def markdown_report(self, md="pdos-report.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding="utf-8") as fout:
            fout.write("# 投影态密度图\n")
            fout.write("**指定能量范围数据图\n")
            fout.write("![pdos-range](./pdos-specified-range.png)\n")
            fout.write("![pdos-atom-range](./pdos-atomproj-specified-range.png)\n")
            fout.write("![tdos-range](./tdos-specified-range.png)\n")
            fout.write("**所有可获取能量范围数据图**\n")
            fout.write("![pdos-all](./pdos-all-energy-available.png)\n")
            fout.write("![pdos-atom-all](./pdos-atomproj-all-energy-available.png)\n")
            fout.write("![tdos-all](./tdos-all-energy-available.png)\n")

    def export(self, directory="tmp-qe-static", plotrange=[0, 1.0], atomtoproj=[], fontsize=10):
        os.chdir(directory)
        os.system("mkdir -p post-processing")
        self.plot_elem_orb_l_proj(plotrange=plotrange, filename="post-processing/pdos-specified-range.png", fontsize=fontsize)
        self.plot_atom_orb_l_proj(plotrange=plotrange, atomtoproj=atomtoproj, filename="post-processing/pdos-atomproj-specified-range.png", fontsize=fontsize)
        self.plot_tdos(plotrange=plotrange, filename="post-processing/tdos-specified-range.png")
        # also plot the all data
        self.plot_elem_orb_l_proj(plotrange=[0, 1.0], filename="post-processing/pdos-all-energy-available.png", fontsize=fontsize)
        self.plot_atom_orb_l_proj(plotrange=[0, 1.0], atomtoproj=atomtoproj, filename="post-processing/pdos-atomproj-all-energy-available.png", fontsize=fontsize)
        self.plot_tdos(plotrange=[0, 1.0], filename="post-processing/tdos-all-energy-available.png")
        self.markdown_report(md="post-processing/pdos-report.md")
        os.chdir("../")
    #



class pdos_post:
    """
    """
    def __init__(self):
        self.data = {} # contain the pdos data but not the tdos
        self.tdos = None
        self.energies = None

    def get_data(self, directory="tmp-qe-static", filpdos="projwfc"):
        """
        this function first try to get fermi energy from the nscfout file
        and if nscfout doesn't exist it will try to extract fermi energy
        from scfout. if both don't exist, it will stop the program and
        print out the warnnig, which guarantee that the fermi energy is
        always shifted to 0

        atomorb is a string in format like this atm#1(Li)_wfc#2(s).
        """
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("pdos post:\n")
            print("  directory of previous scf calculattion not found!\n")
            sys.exit(1)

        os.chdir(directory)
        os.system("ls | grep '%s.pdos_' > projwfc-pdos-file.data" % filpdos)
        with open("projwfc-pdos-file.data", 'r') as fin:
            for line in fin:
                if line.split(".")[1] == "pdos_tot\n":
                    with open(line.split()[0], 'r') as f:
                        f.readline()
                        self.tdos = np.loadtxt(f)
                    continue
                atmorb = line.split("_")[1]+"_"+line.split("_")[2].split()[0]
                # atomorb is a string in format like this atm#1(Li)_wfc#2(s).
                with open(line.split()[0], 'r') as f:
                    f.readline()
                    self.data[atmorb] = np.loadtxt(f)
        os.chdir("../")

        self.energies = self.tdos[:, 0]

        # get fermi energy from nscf output
        scfout = "static-scf.out"
        nscfout = "static-nscf.out"
        if os.path.exists(os.path.join(directory, nscfout)):
            with open(os.path.join(directory, nscfout), 'r') as fin:
                for line in fin:
                    if len(line.split()) == 0:
                        continue
                    if line.split()[0] == "the" and line.split()[1] == "Fermi":
                        efermi = float(line.split()[4])
        elif os.path.exists(os.path.join(directory, scfout)):
            with open(os.path.join(directory, scfout), 'r') as fin:
                for line in fin:
                    if len(line.split()) == 0:
                        continue
                    if line.split()[0] == "the" and line.split()[1] == "Fermi":
                        efermi = float(line.split()[4])
        else:
            print("===========================================================\n")
            print("                Warning !!!\n")
            print("===========================================================\n")
            print("PDOS postprocessing:\n")
            print("must provide nscfout or at least scfout to get Fermi energy\n")
            sys.exit(1)
        # shift fermie energy to 0
        for i in range(len(self.energies)):
            self.energies[i] = self.energies[i] - efermi
        print("===============================================\n")
        print("qe.post.pdos:\n")
        print("we automatically shift the fermi energy\n")
        print("from %f to 0\n" % efermi)
        print("efermi is read from static-nscf.out\n")
        print("or statis-scf.out, if static-nscf.out is not available\n")
        #

    def plot_elem_orb_proj(self, plotrange=[0.0, 1.0], filename="pdos-projected-to-element-and-orbital.png", fontsize=10):
        """
        plotrange:
            a list of two values(between 0 and 1) defining the percentage
            of data to plot.
            plotrange[0]: left boundary of the data to plot
            plotrange[1]: right boundary of the data to plot
            default is plotrange[0] = 0, plotrange[1], in which case
            all the data will be plot.
        """
        data = {}
        for atmorb in self.data:
            key = self.get_elem_type(atmorb)+"-"+self.get_orb_type(atmorb)
            if key in data:
                data[key] = data[key] + self.data[atmorb][:, 2]
            else:
                data[key] = self.data[atmorb][:, 2]

        # plot the pdos in the specified percentage range
        begin = int(len(self.energies)*plotrange[0])
        end = int(len(self.energies)*plotrange[1])
        for key in data:
            plt.plot(self.energies[begin:end], data[key][begin:end], label=key)

        # plot the total dos in the specified percentage range
        plt.plot(self.energies[begin:end], self.tdos[begin:end, 2], label="Total-DOS")

        # set formats
        font = {'size': fontsize}
        plt.tick_params(labelsize=fontsize)

        plt.grid(which="major", axis="x", linewidth=0.75, linestyle="-", color="0.75")
        plt.grid(which="major", axis="y", linewidth=0.75, linestyle="-", color="0.75")
        plt.title("Projected Density of States")
        plt.xlabel(r"$\mathit{E}-\mathit{E}_\mathrm{f} \mathrm{(eV)}$", font)
        plt.ylabel("States", font)
        plt.legend(prop=font)
        plt.tight_layout()
        plt.savefig("%s" % filename)
        plt.close()


    def plot_atom_orb_proj(self, atomtoproj=[], plotrange=[0.0, 1.0], filename="pdos-projected-to-atom-and-orbital.png", fontsize=10):
        """
        plotrange:
            a list of two values(between 0 and 1) defining the percentage
            of data to plot.
            plotrange[0]: left boundary of the data to plot
            plotrange[1]: right boundary of the data to plot
            default is plotrange[0] = 0, plotrange[1], in which case
            all the data will be plot.
        atomtoproj:
            the list of atoms to do the projection. atom number starts with 1
        """
        # plot the data in the specified percentage range
        begin = int(len(self.energies)*plotrange[0])
        end = int(len(self.energies)*plotrange[1])

        # atom projected dos
        for atmorb in self.data:
            if self.get_atom_num(atmorb) in atomtoproj:
                plt.plot(self.energies[begin:end], self.data[atmorb][begin:end, 2], label="Atom(%d):%s-%s" % (self.get_atom_num(atmorb), self.get_elem_type(atmorb), self.get_orb_type(atmorb)))

        # plot the total dos in the specified percentage range
        plt.plot(self.energies[begin:end], self.tdos[begin:end, 2], label="Total-DOS")
        #

        # set formats
        font = {'size': fontsize}
        plt.tick_params(labelsize=fontsize)

        plt.grid(which="major", axis="x", linewidth=0.75, linestyle="-", color="0.75")
        plt.grid(which="major", axis="y", linewidth=0.75, linestyle="-", color="0.75")
        plt.title("Projected(Atom) Density of States")
        plt.xlabel(r"$\mathit{E}-\mathit{E}_\mathrm{f} \mathrm{(eV)}$", font)
        plt.ylabel("States", font)
        plt.legend(prop=font)
        plt.tight_layout()
        plt.savefig("%s" % filename)
        plt.close()


    def plot_tdos(self, plotrange=[0, 1.0], filename="total-dos.png"):
        """
        plotrange:
            a list of two values(between 0 and 1) defining the percentage
            of data to plot.
            plotrange[0]: left boundary of the data to plot
            plotrange[1]: right boundary of the data to plot
            default is plotrange[0] = 0, plotrange[1], in which case
            all the data will be plot.
        """
        # plot the total dos in the specified percentage range
        begin = int(len(self.energies)*plotrange[0])
        end = int(len(self.energies)*plotrange[1])
        #plt.plot(self.energies, self.tdos[:, 2], label="total-dos")
        plt.plot(self.energies[begin:end], self.tdos[begin:end, 2], label="Total-DOS")

        plt.grid(which="major", axis="x", linewidth=0.75, linestyle="-", color="0.75")
        plt.grid(which="major", axis="y", linewidth=0.75, linestyle="-", color="0.75")
        plt.title("Total Density of States")
        plt.xlabel(r"$\mathit{E}-\mathit{E}_\mathrm{f} (eV)$")
        plt.ylabel("States")
        plt.legend()
        plt.tight_layout()
        plt.savefig("%s" % filename)
        plt.close()

    def get_elem_type(self, atmorb):
        """
        get element name from atmorb
        atmorb is the key in self.data
        it's like this:
            atm#1(Li)_wfc#2(s)
        return value of the above input
        will be:
            'Li'
        """
        return atmorb.split("(")[1].split(")")[0]

    def get_orb_type(self, atmorb):
        """
        get element name and orb from atmorb
        atmorb is the key in self.data
        it's like this:
            atm#1(Li)_wfc#2(s)
        return value of the above input
        will be:
            '2(s)'
        """
        return atmorb.split("#")[2]

    def get_atom_num(self, atmorb):
        """
        get atom name from atmorb
        atmorb is the key in self.data
        it's like this:
            atm#1(Li)_wfc#2(s)
        return value of the above input
        will be:
            1
        """
        return int(atmorb.split("(")[0].split("#")[1])

    def markdown_report(self, md="pdos-report.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding="utf-8") as fout:
            fout.write("# 投影态密度图\n")
            fout.write("**指定能量范围数据图\n")
            fout.write("![pdos-range](./pdos-specified-range.png)\n")
            fout.write("![pdos-atom-range](./pdos-atomproj-specified-range.png)\n")
            fout.write("![tdos-range](./tdos-specified-range.png)\n")
            fout.write("**所有可获取能量范围数据图**\n")
            fout.write("![pdos-all](./pdos-all-energy-available.png)\n")
            fout.write("![pdos-atom-all](./pdos-atomproj-all-energy-available.png)\n")
            fout.write("![tdos-all](./tdos-all-energy-available.png)\n")

    def export(self, directory="tmp-qe-static", plotrange=[0, 1.0], atomtoproj=[], fontsize=10):
        os.chdir(directory)
        os.system("mkdir -p post-processing")
        self.plot_elem_orb_proj(plotrange=plotrange, filename="post-processing/pdos-specified-range.png", fontsize=fontsize)
        self.plot_atom_orb_proj(plotrange=plotrange, atomtoproj=atomtoproj, filename="post-processing/pdos-atomproj-specified-range.png", fontsize=fontsize)
        self.plot_tdos(plotrange=plotrange, filename="post-processing/tdos-specified-range.png")
        # also plot the all data
        self.plot_elem_orb_proj(plotrange=[0, 1.0], filename="post-processing/pdos-all-energy-available.png", fontsize=fontsize)
        self.plot_atom_orb_proj(plotrange=[0, 1.0], atomtoproj=atomtoproj, filename="post-processing/pdos-atomproj-all-energy-available.png", fontsize=fontsize)
        self.plot_tdos(plotrange=[0, 1.0], filename="post-processing/tdos-all-energy-available.png")
        self.markdown_report(md="post-processing/pdos-report.md")
        os.chdir("../")
    #
