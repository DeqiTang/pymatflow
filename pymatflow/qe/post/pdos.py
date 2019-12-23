#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

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

    def plot_elem_orb_proj(self, plotrange=[0.0, 1.0], filename="pdos-projected-to-element-and-orbital.png"):
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
        
        plt.grid(which="major", axis="x", linewidth=0.75, linestyle="-", color="0.75")
        plt.grid(which="major", axis="y", linewidth=0.75, linestyle="-", color="0.75")
        plt.title("Projected Density of States")
        plt.xlabel(r"$\mathit{E}-\mathit{E}_\mathrm{f} \mathrm{(eV)}$")
        plt.ylabel("States")
        plt.legend()
        plt.tight_layout()
        plt.savefig("%s" % filename)
        plt.close()


    def plot_atom_orb_proj(self, atomtoproj=[], plotrange=[0.0, 1.0], filename="pdos-projected-to-atom-and-orbital.png"):
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

        plt.grid(which="major", axis="x", linewidth=0.75, linestyle="-", color="0.75")
        plt.grid(which="major", axis="y", linewidth=0.75, linestyle="-", color="0.75")
        plt.title("Projected(Atom) Density of States")
        plt.xlabel(r"$\mathit{E}-\mathit{E}_\mathrm{f} \mathrm{(eV)}$")
        plt.ylabel("States")
        plt.legend()
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

    def export(self, directory="tmp-qe-static", plotrange=[0, 1.0], atomtoproj=[]):
        os.chdir(directory)
        self.plot_elem_orb_proj(plotrange=plotrange, filename="pdos-specified-range.png")
        self.plot_atom_orb_proj(plotrange=plotrange, atomtoproj=atomtoproj, filename="pdos-atomproj-specified-range.png")
        self.plot_tdos(plotrange=plotrange, filename="tdos-specified-range.png")
        # also plot the all data
        self.plot_elem_orb_proj(plotrange=[0, 1.0], filename="pdos-all-energy-available.png")
        self.plot_atom_orb_proj(plotrange=[0, 1.0], atomtoproj=atomtoproj, filename="pdos-atomproj-all-energy-available.png")
        self.plot_tdos(plotrange=[0, 1.0], filename="tdos-all-energy-available.png")
        self.markdown_report()
        os.chdir("../")
    #
