#!/usr/bin/bash
# _*_ coding: utf-8 _*_

import os
import matplotlib.pyplot as plt
import sys

class converge_post:
    def __init__(self):
        # analyse the result
        self.criteria_for_ecutwfc = 7.35e-4 # 10 meV = 7.35e-4 Ry
        self.criteria_for_ecutrho = 7.35e-4
        self.criteria_for_kpoints = 7.35e-3 # 7.35e-4

    def postprocess(self, directory, converge):
        """
        directory:
            directory of the converge test running
        converge:
            type of the converge running, it can be
            ecutwfc, ecutrho, kpoints
        """
        os.chdir(directory)
        if converge == "ecutwfc":
            os.system("ls | grep '.out' | grep 'ecutwfc-' > outfiles.data")
        elif converge == "ecutrho":
            os.system("ls | grep '.out' | grep 'ecutrho-' > outfiles.data")
        elif converge == "kpoints":
            os.system("ls | grep '.out' | grep 'kpoints-' > outfiles.data")
        else:
            print("qe.post.converge.converge_post class can only deal with\n")
            print("ecutwfc, ecutrho, kpoints now\n")
            sys.exit(1)

        outfiles = []
        x_all = []
        with open("outfiles.data", 'r') as fin:
            for line in fin:
                outfiles.append(line.split("\n")[0])
                x_all.append(int(line.split(".")[0].split("-")[1]))
        
        if os.path.exists("energy-x.data"):
            # there exists and old energy-x.data, remove it and generate the new one
            os.system("rm energy-x.data")

        for f in outfiles:
            os.system("cat %s | grep '!    total energy' >> energy-x.data" % f)

        energy_all = []
        with open("energy-x.data", 'r') as fin:
            for line in fin:
                energy_all.append(float(line.split()[4]))

        plt.plot(x_all, energy_all, marker='o')
        if converge == "ecutwfc":
            plt.title("Ecutwfc Converge Test", fontweight='bold', color='red')
            plt.xlabel("Ecutwfc (Ry)")
            plt.ylabel("Energy (Ry)")
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-ecutwfc.png")
            self.md_report(converge=converge, suggested=x_all[self.judge(energy_all, self.criteria_for_ecutwfc)])
        elif converge == "ecutrho":
            plt.title("Ecutrho COnverge Test", fontweight='bold', color='red')
            plt.xlabel("Ecurho (Ry)")
            plt.ylabel("Energy (Ry)")        
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-ecutrho.png")
            self.md_report(converge=converge, suggested=x_all[self.judge(energy_all, self.criteria_for_ecutrho)])
        elif converge == "kpoints":
            plt.title("Kpoints Converge Test", fontweight='bold', color='red')
            plt.xlabel("Kpoints")
            plt.ylabel("Energy (Ry)")
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-kpoints.png")
            self.md_report(converge=converge, suggested=x_all[self.judge(energy_all, self.criteria_for_kpoints)])
        else:
            print("qe.post.converge.converge_post class can only deal with\n")
            print("ecutwfc, ecutrho, kpoints now\n")
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
