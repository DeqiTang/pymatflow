#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import numpy as np
import matplotlib.pyplot as plt

class dos_post:
    """
    Note:
        there are three columns in xxx.dos
        E (eV)   dos(E)     Int dos(E)
        I plot E-dos(E), but don't know what
        Int dos(E) represents(so leave it out
        now)

        automatically shift fermi energy to 0
    """
    def __init__(self):
        self.efermi = None
        self.data = None

    def get_data(self, directory="tmp-qe-static", fildos="dosx.dos"):
        os.chdir(directory)
        with open(fildos, 'r') as fin:
            self.efermi = float(fin.readline().split()[8]) # EFermi in eV
            self.data = np.loadtxt(fin)
        os.chdir("../")
        self.shift_efermi()

    def shift_efermi(self):
        """
        shift efermi to 0
        """
        self.data[:, 0] = self.data[:, 0] - self.efermi

    def plot_dos(self, image="dosx.png"):
        plt.plot(self.data[:, 0], self.data[:, 1], label="TDOS")
        plt.vlines(0, 0, 10, linestyle="dashed", label="Fermi Energy")
        plt.legend()
        plt.savefig("%s" % image)
        plt.show()
    
    def export(self, directory="tmp-qe-static"):
        os.chdir(directory)
        self.plot_dos()
        os.chdir("../")
