#!/usr/bin/bash
# _*_ coding: utf-8 _*_

import os
import matplotlib.pyplot as plt
import sys

class converge_post:
    def __init__(self):
    # analyse the result
        pass

    def ecutwfc(self, directory="tmp-qe-ecutwfc"):
        os.chdir(directory)
        os.system("ls | grep '.out' | grep 'ecutwfc-'> outfiles.data")
        outfiles = []
        ecutwfc_all = []
        with open("outfiles.data", 'r') as fin:
            for line in fin:
                outfiles.append(line.split("\n")[0])
                ecutwfc_all.append(int(line.split(".")[0].split("-")[1]))
        
        
        for f in outfiles:
            os.system("cat %s | grep '!    total energy' >> energy-ecutwfc.data" % f)

        energy_all = []
        with open("energy-ecutwfc.data", 'r') as fin:
            for line in fin:
                energy_all.append(float(line.split()[4]))

        plt.plot(ecutwfc_all, energy_all, marker='o')
        plt.title("Ecutwfc Converge Test", fontweight='bold', color='red')
        plt.xlabel("Ecutwfc (Ry)")
        plt.ylabel("Energy (Ry)")
        plt.tight_layout()
        plt.grid(True)
        plt.savefig("energy-ecutwfc.png")
        plt.show()
        os.chdir("../")


