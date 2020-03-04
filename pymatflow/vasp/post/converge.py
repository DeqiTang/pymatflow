#!/usr/bin/bash
# _*_ coding: utf-8 _*_

import os
import re
import matplotlib.pyplot as plt
import sys

class converge_post:
    def __init__(self):
        # analyse the result
        self.criteria_for_encut = 0.01 # 10 meV
        self.criteria_for_sigma = 0.01
        self.criteria_for_kpoints = 0.01

    def postprocess(self, directory, converge):
        """
        directory:
            directory of the converge test running
        converge:
            type of the converge running, it can be
            encut, sigma, kpoints
        """
        os.chdir(directory)
        if converge == "encut":
            x_all = []
            for f in os.listdir():
                if f.split("-")[0] == "encut" and len(f.split("-")) == 2:
                    x_all.append(int(f.split("-")[1]))
            # we must sort the x_all
            x_all.sort()

            outfiles = []
            for x in x_all:
                outfiles.append("encut-%d/OUTCAR" % x)

        elif converge == "sigma":
            x_all = []
            for f in os.listdir():
                if f.split("-")[0] == "sigma" and len(f.split("-")) == 2:
                    x_all.append(float(f.split("-")[1]))
            # we must sort the x_all
            x_all.sort()

            outfiles = []
            for x in x_all:
                outfiles.append("sigma-%.6f/OUTCAR" % x)

        elif converge == "kpoints":
            x_all = []
            for f in os.listdir():
                if f.split("-")[0] == "kpoints" and len(f.split("-")) == 2:
                    x_all.append(int(f.split("-")[1]))
            # we must sort the x_all
            x_all.sort()

            outfiles = []
            for x in x_all:
                outfiles.append("kpoints-%d/OUTCAR" % x)
        else:
            print("vasp.post.converge.converge_post class can only deal with\n")
            print("encut, sigma, kpoints now\n")
            sys.exit(1)

        
        if os.path.exists("energy-x.data"):
            # there exists and old energy-x.data, remove it and generate the new one
            os.system("rm energy-x.data")

        for f in outfiles:
            os.system("cat %s | grep 'energy  without entropy=' >> energy-x.data" % f)

        energy_all = []
        with open("energy-x.data", 'r') as fin:
            for line in fin:
                energy_all.append(float(line.split()[3]))

        plt.plot(x_all, energy_all, marker='o')
        if converge == "encut":
            plt.title("Encut Converge Test", fontweight='bold', color='red')
            plt.xlabel("Encut (eV)")
            plt.ylabel("Energy (eV)")
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-encut.png")
            self.md_report(converge=converge, suggested=x_all[self.judge(energy_all, self.criteria_for_encut)])
        elif converge == "sigma":
            plt.title("Sigma Converge Test", fontweight='bold', color='red')
            plt.xlabel("Sigma (eV)")
            plt.ylabel("Energy (eV)")        
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-sigma.png")
            self.md_report(converge=converge, suggested=x_all[self.judge(energy_all, self.criteria_for_sigma)])
        elif converge == "kpoints":
            plt.title("Kpoints Converge Test", fontweight='bold', color='red')
            plt.xlabel("Kpoints")
            plt.ylabel("Energy (eV)")
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-kpoints.png")
            self.md_report(converge=converge, suggested=x_all[self.judge(energy_all, self.criteria_for_kpoints)])
        else:
            print("vasp.post.converge.converge_post class can only deal with\n")
            print("encut, sigma, kpoints now\n")
            sys.exit(1)
        
        plt.show()
        os.chdir("../")

    #

    def judge(self, energies, criteria):
        #
        # return the index of the recommended value
        # if criteria not reached return -1 which
        # means the last params
        for i in range(1, len(energies)):
            deltae = energies[i] - energies[i-1]
            if abs(deltae) < criteria:
                return i
        # didn't satisfy the criteria return the last index
        return -1

    def md_report(self, converge, suggested):
        with open("%s-test-report.md" % converge, 'w') as fout:
            fout.write("## %s 收敛测试\n" % converge)
            fout.write("推荐值:\n")
            fout.write("```\n")
            fout.write("%s: %s\n" % (converge, str(suggested)))
            fout.write("```\n")
            fout.write("注意:\n")
            fout.write("- 这个推荐值是你的测试范围内的推荐值\n")
            fout.write("- 如果该值是测试范围的最后一个值, 可能意味着没有达到收敛判剧\n")
            if converge == "encut":
                fout.write("- 推荐值的依据为%s前后两值能量差小于%f[eV]\n" % (converge, self.criteria_for_encut))
            if converge == "sigma":
                fout.write("- 推荐值的依据为%s前后两值能量差小于%f[eV]\n" % (converge, self.criteria_for_sigma))
            if converge == "kpoints":
                fout.write("- 推荐值的依据为%s前后两值能量差小于%f[eV]\n" % (converge, self.criteria_for_kpoints))
            fout.write("### 能量变化趋势图\n")
            fout.write("![energy-x](energy-%s.png)\n" % converge)
