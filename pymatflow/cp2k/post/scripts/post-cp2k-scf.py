#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import datetime
import argparse
import matplotlib.pyplot as plt

from pymatflow.cp2k.post.scf import scf_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of static SCF running", type=str, default="tmp-cp2k-static")
    parser.add_argument("-f", "--file", help="output of scf running: name only without pasth", type=str, default="static-scf.out")

    args = parser.parse_args()

    os.chdir(args.directory)
    os.system("mkdir -p post-processing")
    os.chdir("post-processing")
    scf = scf_out()
    scf.get_info(file="../%s" % args.file)

    # now output the scf information in the current directory(post-processing

    plt.plot(scf.run_info["scf_energies"])
    plt.title("Total energy per scf step")
    plt.xlabel("Scf step")
    plt.ylabel("Total energy")
    plt.tight_layout()
    plt.savefig("total-energy-per-scf-step.png")
    plt.close()
    pass

    with open("scf-info.md", 'w', encoding='utf-8') as fout:
        fout.write("# SCF实验统计\n")
        fout.write("## SCF参数\n")
        for item in scf.scf_params:
            fout.write("- %s: %s\n" % (item, str(scf.scf_params[item])))
        fout.write("## 运行信息\n")
        # calculate the running time and print it out
        # note that the accuracy for the seconds is not fully guranteed
        # e.g. 2019-11-26 12:09:36.487 is read as 2019-11-26 12:09:36
        start = datetime.datetime.strptime(scf.run_info["start_time"].split()[7]+"-"+scf.run_info["start_time"].split()[8].split(".")[0], "%Y-%m-%d-%H:%M:%S")
        stop = datetime.datetime.strptime(scf.run_info["stop_time"].split()[7]+"-"+scf.run_info["stop_time"].split()[8].split(".")[0], "%Y-%m-%d-%H:%M:%S")
        delta_t = stop -start
        fout.write("- Time consuming:\n")
        fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
        # end the time information
        for item in scf.run_info:
            fout.write("- %s: %s\n" % (item, str(scf.run_info[item])))

        fout.write("## 运行信息图示\n")
        fout.write("![dft-energy-per-scf](total-energy-per-scf-step.png)")
    os.chdir("../")
    os.chdir("../")


    # --------------------------------------------------------------------------
    # print information to the terminal
    # --------------------------------------------------------------------------
    print("=====================================================================\n")
    print("                           post-cp2k-scf.py\n")
    print("---------------------------------------------------------------------\n")
    print("\n")
