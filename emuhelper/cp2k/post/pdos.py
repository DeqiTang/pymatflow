#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
from scipy.ndimage.filters import gaussian_filter1d as gf
from scipy import interpolate

import matplotlib.pyplot as plt

class element_pdos:
    """
    """
    def __init__(self, pdos_file):
        self.pdos_file = pdos_file
        self.kind = None
        self.fermi = None
        self.orbitals = {}
        self.data = None
        self.get_data()

    def get_data(self):
        with open(self.pdos_file, 'r') as fin:
            line = fin.readline()
            self.kind = line.split()[6]
            self.fermi = float(line.split()[15]) # a.u.
            line = fin.readline()
            for item in line.split()[5:]:
                self.orbitals[item] = None
            self.data = np.loadtxt(fin)
        #
        j = 3
        for orb in self.orbitals:
            self.orbitals[orb] = self.data[:, j]
            j += 1
        

class pdos_post:
    """
    an abstraction of the .PDOS file from cp2k running
    """
    def __init__(self, pdos_files, step=0.001, sigma=0.01):
        self.data = []
        for f in pdos_files:
            self.data.append(element_pdos(f))

        self.step = step
        self.sigma = sigma

    def convolute(self, x, y, step=0.001, sigma=0.03):
        """
        need to do a convolution plot using gaussian to get a smooth DOS.
        Reference:
            https://matthew-brett.github.io/teaching/smoothing_as_convolution.html
        maybe we can use scipy.ndimage.gaussian_filter1d to do it conveniently
        """
        nx = int((max(x) - min(x)) / step) + 1
        xnew = np.linspace(min(x), max(x), nx)
        ynew = np.zeros(nx)
        for i in range(nx):
            xi = min(x) + i * step # xi is the center of the i gaussian kernel
            kernel_at_i = np.exp(-(x - xi) ** 2 / (2 * sigma ** 2)) / (sigma*np.sqrt(2.0*np.pi))
            ynew[i] = y.dot(kernel_at_i)
        
        return xnew, ynew

    def plot_tdos(self):
        energy = self.data[0].data[:, 1]
        tdos = np.zeros(len(energy))
        for i in self.data:
            for j in range(3, i.data.shape[1]):
                tdos += i.data[:, j]
        
        x, y = self.convolute(energy, tdos, step=self.step, sigma=self.sigma)
        plt.plot(x, y)
        plt.title("Total DOS(conv)")
        plt.savefig("tdos.png")
        plt.close()

   
    def export(self):
        self.plot_tdos()

