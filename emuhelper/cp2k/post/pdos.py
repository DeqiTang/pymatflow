#!/usr/bin/env python
# _*_ coding: utf-8 _*_

class pdos:
    """
    an abstraction of the .PDOS file from cp2k running
    """
    def __init__(self, fname):
        self.eigenvalue = []
        self.occupation = []
        self.orbitals = {}
        with open(fname, 'r') as fin:
            line = fin.readline()
            self.element = line.split()[6]
            self.iterstep = int(line.split()[12].split(",")[0])
            self.fermi = float(line.split()[15]) # a.u.
            line = fin.readline()
            for item in line.split()[5:]:
                self.orbitals[item] = []
            for line in fin:
                self.eigenvalue.append(float(line.split()[1]))
                self.occupation.append(float(line.split()[2]))
                i = 3
                for item in self.orbitals:
                    self.orbitals[item].append(float(line.split()[i]))
                    i += 1
            # end file
        self.energy = [(i - self.fermi) * 27.211384523 for i in self.eigenvalue]
        self.tpdos = []
        for i in range(len(self.energy)):
            self.tpdos.append(float(0))
            for item in self.orbitals:
                self.tpdos[i] = self.tpdos[i] + float(self.orbitals[item][i])

    def help(self):
        pass
