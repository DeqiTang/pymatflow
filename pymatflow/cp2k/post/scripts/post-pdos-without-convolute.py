#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import matplotlib.pyplot as plt

fermi = None
energies = []
orbitals = {}
with open(sys.argv[1], 'r') as fin:
    line1 = fin.readline()
    fermi = line1.split()[15]
    line2 = fin.readline()
    for item in line2.split()[5:]:
        orbitals[item] = []

    n_orbitals = len((orbitals))
    for line in fin:
        energies.append(float(line.split()[1]))
        i = 3
        for orbital in orbitals:
            orbitals[orbital].append(line.split()[i])
            i += 1

all_orbitals = []
for i in range(len(energies)):
    all_orbitals.append(float(0))
    for item in orbitals:
        all_orbitals[i] = all_orbitals[i] + float(orbitals[item][i])
#
plt.plot(energies, all_orbitals)
plt.show()
