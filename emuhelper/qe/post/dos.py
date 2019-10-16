#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
import matplotlib.pyplot as plt

class qe_dos:
    """
    Note:
        there are three columns in xxx.dos
        E (eV)   dos(E)     Int dos(E)
        I plot E-dos(E), but don't know what
        Int dos(E) represents(so leave it out
        now)

        automatically shift fermi energy to 0
    """
    def __init__(self, dosfile):
        self.dosfile = dosfile
        self.efermi = None
        self.data = None
        self.load_data()
        self.shift_efermi()

    def load_data(self):
        with open(self.dosfile, 'r') as fin:
            self.efermi = float(fin.readline().split()[8]) # EFermi in eV
            self.data = np.loadtxt(fin)

    def shift_efermi(self):
        """
        shift efermi to 0
        """
        self.data[:, 0] = self.data[:, 0] - self.efermi

    def plot_dos(self):
        plt.plot(self.data[:, 0], self.data[:, 1], label="TDOS")
        plt.vlines(0, 0, 10, linestyle="dashed", label="Fermi Energy")
        plt.legend()
        plt.savefig("%s.png" % self.dosfile)
        plt.show()

