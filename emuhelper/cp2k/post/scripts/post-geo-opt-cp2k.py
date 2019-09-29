#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import matplotlib.pyplot as plt

"""
usage:
    post-geo-opt-cp2k.py xxx
    xxx is the output file of the GEO_OPT run
"""



os.system("cat %s | grep 'ENERGY| Total FORCE_EVAL' > energy-per-geo-step.data" % (sys.argv[1]))

os.system("cat %s | grep '*** SCF run converged in' > scf-steps.data" % (sys.argv[1]))

energies = []
with open("energy-per-geo-step.data", 'r') as fin:
    for line in fin:
        energies.append(float(line.split()[8]))

scf_steps = []
with open("scf-steps.data", 'r') as fin:
    for line in fin:
        scf_steps.append(int(line.split()[5]))

ion_steps = [i for i in range(len(energies))]
plt.plot(ion_steps, energies)
plt.show()
