#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import argparse
import copy
import datetime
import subprocess
import matplotlib.pyplot as plt

from pymatflow.abinit.post.scf import scf_out



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously static running directory", type=str, default="tmp-abinit-static")
    parser.add_argument("--scfout", help="output file of static calculation", type=str, default="static-scf.out")

    args = parser.parse_args()

    #os.chdir(args.directory)
    #task = scf_post(args.scfout)
    #task.export()
    #os.chdir("../")



    os.chdir(args.directory)
    scf = scf_out()
    scf.get_info(file=args.scfout)

    os.system("mkdir -p post-processing")
    os.chdir("post-processing")

    #plt.plot(self.run_info["iterations"])
    #plt.title("Iterations per SCF")
    #plt.xlabel("Scf cycles ")
    #plt.ylabel("iterations")
    #plt.tight_layout()
    #plt.savefig("iterations-per-scf.png")
    #plt.close()


    with open("scf-info.md", 'w', encoding='utf-8') as fout:
        fout.write("# 静态计算实验统计\n")
        fout.write("## 计算参数\n")
        for item in self.scf_params:
            fout.write("- %s: %s\n" % (item, str(self.scf_params[item])))
        fout.write("## 运行信息\n")
        # calculate the running time and print it out

        # end the time information
        for item in self.run_info:
            fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

        fout.write("## 运行信息图示\n")

        fout.write("Iterations per SCF\n")
        fout.write("![Iterations per SCF](iterations-per-scf.png)\n")

    os.chdir("../")
    os.chdir("../")


    # --------------------------------------------------------------------------
    # print information to the terminal
    # --------------------------------------------------------------------------
    print("=====================================================================\n")
    print("                           post-abinit-scf.py\n")
    print("---------------------------------------------------------------------\n")
    print("\n")
