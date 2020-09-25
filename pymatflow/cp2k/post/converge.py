#!/usr/bin/bash
# _*_ coding: utf-8 _*_

import os
import re
import matplotlib.pyplot as plt
import sys

class converge_post:
    def __init__(self):
        # analyse the result
        self.criteria_for_cutoff = 3.67e-4 # in unit of Hartree, 10 meV = 0.00036749308136648884 Hartree
        self.criteria_for_rel_cutoff = 3.67e-4
        self.criteria_for_kpoints = 3.67e-4


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
        with open("%s-test-report.md" % converge.lower(), 'w') as fout:
            fout.write("## %s 收敛测试\n" % converge.upper())
            fout.write("推荐值:\n")
            fout.write("```\n")
            fout.write("%s: %s\n" % (converge.upper(), str(suggested)))
            fout.write("```\n")
            fout.write("注意:\n")
            fout.write("- 这个推荐值是你的测试范围内的推荐值\n")
            fout.write("- 如果该值是测试范围的最后一个值, 可能意味着没有达到收敛判剧\n")
            if converge.upper() == "CUTOFF":
                fout.write("- 推荐值的依据为%s前后两值能量差小于%.9f[Hartree]=%.9f[eV]\n" % (converge.upper(), self.criteria_for_cutoff, self.criteria_for_cutoff*27.211396641308))
            if converge.upper() == "REL_CUTOFF":
                fout.write("- 推荐值的依据为%s前后两值能量差小于%.9f[Hartree]=%.9f[eV]\n" % (converge.upper(), self.criteria_for_rel_cutoff, self.criteria_for_rel_cutoff*27.211396641308))
            fout.write("### 能量变化趋势图\n")
            fout.write("![energy-x](energy-%s.png)\n" % converge.lower())


    def postprocess(self, directory, converge):
        """
        directory:
            directory of the converge test running
        converge:
            type of the converge running, it can be
            CUTOFF(cutoff), REL_CUTOFF(rel_cutoff)
            KPOINTS-MANUAL, KPOINTS-AUTO
            KPOINTS-MANUAL and KPOINTS-AUTO are two different
            ways of kpoints testing. with KPOINTS-MAUAL user
            will specify kpoints manually during converge test,
            while with KPOINTS_AUTO user only need to specify
            the test range by kmin, kmax, step, and converge
            manager will generate the kpoints for test correspondingly
        Warning:
            postprocess for KPOINTS-AUTO works normally now,
            but for KPOINTS-MANUAL it may not work right. so
            be cautious when using it.
        """
        os.chdir(directory)
        if converge.upper() == "CUTOFF":
            x_all = []
            for f in os.listdir():
                if f.split(".")[-1] == "out" and f[0:7] == "cutoff-":
                    x_all.append(int(f.split(".")[0].split("-")[-1]))
            # we must sort the x_all
            x_all.sort()

            outfiles = []
            for x in x_all:
                outfiles.append("cutoff-%d.out" % x)

        elif converge.upper() == "REL_CUTOFF":
            x_all = []
            for f in os.listdir():
                if f.split(".")[-1] == "out" and f[0:11] == "rel-cutoff-":
                    x_all.append(int(f.split(".")[0].split("-")[-1]))
            # we must sort the x_all
            x_all.sort()

            outfiles = []
            for x in x_all:
                outfiles.append("rel-cutoff-%d.out" % x)        

        elif converge.upper() == "KPOINTS_AUTO":
            x_all = []
            for f in os.listdir():
                if f.split(".")[-1] == "out" and f[0:8] == "kpoints-":
                    x_all.append(int(f.split(".")[0].split("-")[-1]))
            # we must sort the x_all
            x_all.sort()

            outfiles = []
            for x in x_all:
                outfiles.append("kpoints-%d.out" % x)

        elif converge.upper() == "KPOINTS_MANUAL":
            x_all = []
            for f in os.listdir():
                if f.split(".")[-1] == "out" and f[0:8] == "kpoints-":
                    x_all.append(f.split(".")[0].split("-")[-1])
            outfiles = []
            for x in x_all:
                outfiles.append("kpoints-%s.out" % (x))

        else:
            print("cp2k.post.converge.converge_post class can only deal with\n")
            print("CUTOFF, REL_CUTOFF, KPONITS_AUTO KPOINTS_MANUAL now\n")
            sys.exit(1)

        # 
        os.system("mkdir -p post-processing")
        os.chdir("post-processing")
        if os.path.exists("energy-x.data"):
            # there exists and old energy-x.data, remove it and generate the new one
            os.system("rm energy-x.data")

        for f in outfiles:
            os.system("cat %s | grep 'Total energy:' >> energy-x.data" % os.path.join("../", f))

        energy_all = []
        with open("energy-x.data", 'r') as fin:
            for line in fin:
                energy_all.append(float(line.split()[2]))
        if converge.upper() == "CUTOFF":
            plt.plot(x_all, energy_all, marker='o')
            plt.title("CUTOFF Converge Test", fontweight='bold', color='red')
            plt.xlabel("CUTOFF (Ry)")
            plt.ylabel("Energy (Ha)")
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-%s.png" % converge.lower())
            self.md_report(converge=converge, suggested=x_all[self.judge(energy_all, self.criteria_for_cutoff)])
        elif converge.upper() == "REL_CUTOFF":
            plt.plot(x_all, energy_all, marker='o')
            plt.title("REL_CUTOFF Converge Test", fontweight='bold', color='red')
            plt.xlabel("REL_CUTOFF (Ry)")
            plt.ylabel("Energy (Ha)")        
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-%s.png" % converge.lower())
            self.md_report(converge=converge, suggested=x_all[self.judge(energy_all, self.criteria_for_rel_cutoff)])
        elif converge.upper() == "KPOINTS-AUTO":
            plt.plot(x_all, energy_all, marker='o')
            #plt.xticks(x_all, ["%dx%dx%d" % (x_all[i], x_all[i], x_all[i]) for i in range(len(x_all))])
            plt.title("Kpoints-auto Converge Test", fontweight='bold', color='red')
            plt.xlabel("Kpoints")
            plt.ylabel("Energy (Ha)")        
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-%s.png" % converge.lower())
            self.md_report(converge=converge, suggested=x_all[self.judge(energy_all, self.criteria_for_kpoints)])
        elif converge.upper() == "KPOINTS-MANUAL":
            plt.plot(energy_all, marker='o')
            plt.xticks(range(len(energy_all)), x_all)
            plt.title("Kpoints-manual Converge Test", fontweight='bold', color='red')
            plt.xlabel("Kpoints")
            plt.ylabel("Energy (Ha)")        
            plt.tight_layout()
            plt.grid(True)
            plt.savefig("energy-%s.png" % converge.lower())
            #self.md_report(converge=converge, suggested=x_all[self.judge(energy_all, self.criteria_for_kpoints)])
        else:
            print("qe.post.converge.converge_post class can only deal with\n")
            print("CUTOFF, REL_CUTOFF, KPOINTS now\n")
            sys.exit(1)
        
        plt.show()
        os.chdir("../../")

    #
