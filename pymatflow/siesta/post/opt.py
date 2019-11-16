#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import matplotlib.pyplot as plt


class opt_post:
    """
    """
    def __init__(self):
        pass

    def analysis(self, directory="tmp-siesta-opt", output="geometric-optimization.out"):
        """
        output:
            output file of the geometric running
        """
        # analyse the results

        os.chdir(directory)
        os.system("cat %s | grep 'siesta: E_KS(eV) =' > energy-per-ion-step.data" % (output))

        energies = []
        with open("energy-per-ion-step.data", 'r') as fin:
            for line in fin:
                energies.append(float(line.split()[3]))

        steps = [i for i in range(len(energies))]
        plt.plot(steps, energies)
        plt.show()
        os.chdir("../")

