#!/usr/bin/bash
# _*_ coding: utf-8 _*_

import os
import re
import matplotlib.pyplot as plt
import sys

class ConvergePost:
    def __init__(self):
        # analyse the result
        self.criteria_for_ecutwfc = 7.35e-4 # 10 meV = 7.35e-4 Ry
        self.criteria_for_ecutrho = 7.35e-5
        self.criteria_for_degauss = 7.35e-5 
        self.criteria_for_kpoints = 7.35e-3 # 7.35e-4

    def postprocess(self, directory, converge):
        """
        directory:
            directory of the converge test running
        converge:
            type of the converge running, it can be
            ecutwfc, ecutrho, degauss, kpoints
        """
        os.chdir(directory)
        if converge == "ecutwfc":
            x_all = []
            for f in os.listdir():
                if f.split(".")[-1] == "out" and f[0:8] == "ecutwfc-":
                    x_all.append(int(f.split(".")[0].split("-")[1]))
            # we must sort the x_all
            x_all.sort()

            outfiles = []
            for x in x_all:
                outfiles.append("ecutwfc-%d.out" % x)

        elif converge == "ecutrho":
            x_all = []
            for f in os.listdir():
                if f.split(".")[-1] == "out" and f[0:8] == "ecutrho-":
                    x_all.append(int(f.split(".")[0].split("-")[1]))
            # we must sort the x_all
            x_all.sort()

            outfiles = []
            for x in x_all:
                outfiles.append("ecutrho-%d.out" % x)

        elif converge == "degauss":
            x_all = []
            for f in os.listdir():
                if f.split(".")[-1] == "out" and f[0:8] == "degauss-":
                    x_all.append(float(f.split(".out")[0].split("-")[1]))
            # we must sort the x_all
            x_all.sort()

            outfiles = []
            for x in x_all:
                outfiles.append("degauss-%f.out" % x)

        elif converge == "kpoints":
            x_all = []
            for f in os.listdir():
                if f.split(".")[-1] == "out" and f[0:8] == "kpoints-":
                    x_all.append(int(f.split(".")[0].split("-")[1]))
            # we must sort the x_all
            x_all.sort()

            outfiles = []
            for x in x_all:
                outfiles.append("kpoints-%d.out" % x)
        else:
            print("qe.post.converge.converge_post class can only deal with\n")
            print("ecutwfc, ecutrho, degauss, kpoints now\n")
            sys.exit(1)

        
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

        os.system("mkdir -p post-processing")
        os.chdir("post-processing")
        if converge == "ecutwfc":
            plt.title("Ecutwfc Converge Test", fontweight='bold', color='red')
            plt.xlabel("Ecutwfc (Ry)")
            plt.ylabel("Energy (Ry)")
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-ecutwfc.png")
            self.md_report(converge=converge, suggested=x_all[self.judge(energy_all, self.criteria_for_ecutwfc)])
        elif converge == "ecutrho":
            plt.title("Ecutrho Converge Test", fontweight='bold', color='red')
            plt.xlabel("Ecurho (Ry)")
            plt.ylabel("Energy (Ry)")        
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-ecutrho.png")
            self.md_report(converge=converge, suggested=x_all[self.judge(energy_all, self.criteria_for_ecutrho)])
        elif converge == "degauss":
            plt.title("Degauss Converge Test", fontweight='bold', color='red')
            plt.xlabel("Degauss (Ry)")
            plt.ylabel("Energy (Ry)")        
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-degauss.png")
            self.md_report(converge=converge, suggested=x_all[self.judge_degauss(energy_all, self.criteria_for_degauss)])        
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
        os.chdir("../../")

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

    def judge_degauss(self, energies, criteria):
        #
        # return the index of the recommended value
        # if criteria not reached return -1 which
        # means the last params
        # 
        # degauss is different from ecutwfc and ecutrho as usually 
        # a small value is better for degauss, while a large value is better for ecutwfc and ecutrho.
        # so we define a unique judge function for degauss.
        for i in range(len(energies)-2, -1, -1):
            deltae = energies[i] - energies[i+1]
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
            if converge == "ecutwfc":
                fout.write("- 推荐值的依据为%s前后两值能量差小于%f[Ry]=%f[eV]\n" % (converge, self.criteria_for_ecutwfc, self.criteria_for_ecutwfc*13.6056923))
            if converge == "ecutrho":
                fout.write("- 推荐值的依据为%s前后两值能量差小于%f[Ry]=%f[eV]\n" % (converge, self.criteria_for_ecutrho, self.criteria_for_ecutrho*13.6056923))
            if converge == "degauss":
                fout.write("- 推荐值的依据为%s前后两值能量差小于%f[Ry]=%f[eV]\n" % (converge, self.criteria_for_degauss, self.criteria_for_degauss*13.6056923))            
            if converge == "kpoints":
                fout.write("- 推荐值的依据为%s前后两值能量差小于%f[Ry]=%f[eV]\n" % (converge, self.criteria_for_kpoints, self.criteria_for_kpoints*13.6056923))
            fout.write("### 能量变化趋势图\n")
            fout.write("![energy-x](energy-%s.png)\n" % converge)
