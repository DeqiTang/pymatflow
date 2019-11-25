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

    def plot_elem_orb_proj(self):
        data = {}
        for atmorb in self.data:
            key = self.get_elem_type(atmorb)+"-"+self.get_orb_type(atmorb)
            if key in data:
                data[key] = data[key] + self.data[atmorb][:, 2]
            else:
                data[key] = self.data[atmorb][:, 2]

        for key in data:
            plt.plot(self.energies, data[key], label=key)
        # plot fermi energy
        # plt.vlines(0, 0, 10, label="Fermi energy")
        # abandoned the above plot but use grid to be more tidy
        plt.grid(which="major", axis="x", linewidth=0.75, linestyle="-", color="0.75")
        plt.grid(which="major", axis="y", linewidth=0.75, linestyle="-", color="0.75")
        plt.title("Projected Density of States")
        plt.xlabel("Energy (eV)")
        plt.ylabel("States")
        plt.legend()
        plt.tight_layout()
        plt.savefig("pdos-projected-to-element-and-orbital.png")
        plt.close()
 
    def plot_tdos(self):
        plt.plot(self.energies, self.tdos[:, 2], label="total-dos")
        # plot fermi energy
        # plt.vlines(0, 0, 10, label="Fermi energy")
        # abandoned the above plot but use grid to be more tidy
        plt.grid(which="major", axis="x", linewidth=0.75, linestyle="-", color="0.75")
        plt.grid(which="major", axis="y", linewidth=0.75, linestyle="-", color="0.75")
        plt.title("Total Density of States")
        plt.xlabel("Energy (eV)")
        plt.ylabel("States")
        plt.legend()
        plt.tight_layout()
        plt.savefig("total-dos.png")
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

    def markdown_report(self, md="pdos-report.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding="utf-8") as fout:
            fout.write("# 投影态密度图\n")
            fout.write("")

    def export(self, directory="tmp-qe-static"):
        os.chdir(directory)
        self.plot_elem_orb_proj()
        self.plot_tdos()
        self.markdown_report()
        os.chdir("../")
    #
