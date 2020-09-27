#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import copy
import numpy as np

import matplotlib.pyplot as plt

class element_pdos:
    """
    """
    def __init__(self, pdos_file):
        self.pdos_file = pdos_file
        self.kind = None
        self.fermi = None # in a.u.
        self.data = None
        self.orbitals = []
        self.get_data()

    def get_data(self):
        with open(self.pdos_file, 'r') as fin:
            line = fin.readline()
            self.kind = line.split()[6]
            self.fermi = float(line.split()[15]) # a.u.
            line = fin.readline()
            for item in line.split()[5:]:
                self.orbitals.append(item)
            self.data = np.loadtxt(fin)
       
        # shift the Fermi energy to 0 and onvert from a.u. to eV
        self.data[:, 1] = self.data[:, 1] - self.fermi
        self.data[:, 1] = self.data[:, 1] * 27.211324570273

        

class pdos_post:
    """
    an abstraction of the .PDOS file from cp2k running
    Note:
        the order of execution convolution and PDOS addition will not affect the result.
        which means we can first do the convolution for every orbital of every element
        and the total dos can be the summation of the convoluted results.
    """
    def __init__(self):
        self.data = []

    def get_data(self, pdos_files):
        for f in pdos_files:
            self.data.append(element_pdos(f))

    def delta(self,emin,emax,npts,energy,width):
        """Return a delta-function centered at energy
        
        Parameters
        ----------
        emin: float
            minimun eigenvalue
        emax: float
            maximun eigenvalue
        npts: int
            Number of points in the smeared pdos
        energy: float
            energy where the gaussian is centered
        width: float
            dispersion parameter

        Return 
        ------
        delta: numpy array
            array of delta function values

        """
        
        energies = np.linspace(emin, emax, npts)
        x = -((energies - energy) / width)**2
        return np.exp(x) / (np.sqrt(np.pi) * width)

    def convolute(self, width=0.2):
        """
        need to do a convolution plot using gaussian to get a smooth DOS.
        self.data_smearing are in the same format as self.data, but it is smoothed by smearing
        """
        self.data_smearing = copy.deepcopy(self.data)
        for i in range(len(self.data)):
            emin = np.min(self.data[i].data[:, 1])
            emax = np.max(self.data[i].data[:, 1])
            for j in range(3, self.data[i].data.shape[1]):
                self.data_smearing[i].data[:,j] = np.zeros(len(self.data_smearing[i].data[:,j]))
                for e, pd in zip(self.data[i].data[:, 1], self.data[i].data[:, j]):
                    self.data_smearing[i].data[:, j] += pd * self.delta(emin, emax, len(self.data[i].data[:, 1]), e, width)
            self.data_smearing[i].data[:, 1] = np.linspace(emin, emax, len(self.data_smearing[i].data[:, 1]))
            
    def export_smearing(self, directory):
        for i in range(len(self.data_smearing)):
            with open(os.path.join(directory, "%s.smearing.data" % os.path.basename(self.data_smearing[i].pdos_file)), 'w') as fout:
                fout.write("#MO\tEigenvalue[eV]\tOccupation")
                for orb in self.data_smearing[i].orbitals:
                    fout.write("\t%s" % orb)
                fout.write("\ttotal(%s)" % self.data_smearing[i].kind) # added total dos for the element
                fout.write("\n")
                fout.write("#int\teV\tfloat")
                for orb in self.data_smearing[i].orbitals:
                    fout.write("\tfloat")
                fout.write("\tfloat") # added total dos for the element
                fout.write("\n")
                fout.write("#Projected DOS for atomic kind %s, eigenvalue are transfered from a.u. to eV already and efermi shifted to 0\n" % self.data_smearing[i].kind)

                for j in range(self.data_smearing[i].data.shape[0]):
                    fout.write("%d" % self.data_smearing[i].data[j, 0])
                    fout.write("\t%.6f\t%.6f" % (self.data_smearing[i].data[j, 1], self.data_smearing[i].data[j, 2]))
                    for k in range(3, self.data_smearing[i].data.shape[1]):
                        fout.write("\t%.8f" % self.data_smearing[i].data[j, k])
                    fout.write("\t%.8f" % np.sum(self.data_smearing[i].data[j, 3:-1]))
                    fout.write("\n")

        with open(os.path.join(directory, "element-projected.pdos.smearing.data"), "w") as fout:
            fout.write("#MO\tEigenvalue[eV]")
            for i in range(len(self.data_smearing)):
                fout.write("\t%s" % self.data_smearing[i].kind)
            fout.write("\tTDOS\n")
            fout.write("#int\teV")
            for i in range(len(self.data_smearing)):
                fout.write("\tfloat")
            fout.write("\tfloat\n")
            fout.write("#eigenvalue are transfered from a.u. to eV already and efermi shifted to 0\n")
            for j in range(self.data_smearing[0].data.shape[0]):
                fout.write("%d" % (j+1))
                fout.write("\t%.6f" % self.data_smearing[0].data[j,1])
                for i in range(len(self.data_smearing)):
                    fout.write("\t%.8f" % np.sum(self.data_smearing[i].data[j, 3:-1]))
                fout.write("\t%.8f\n" % np.sum([np.sum(self.data_smearing[ele].data[j, 3:-1]) for ele in range(len(self.data_smearing))]))
            #