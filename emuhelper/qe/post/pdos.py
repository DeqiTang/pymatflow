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
        plt.title("Projected Density of States")
        plt.xlabel("Energy (eV)")
        plt.ylabel("States")
        plt.legend()
        plt.tight_layout()
        plt.savefig("pdos-projected-to-element-and-orbital.png")
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
        self.markdown_report()
        os.chdir("../")
    #
