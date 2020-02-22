#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import datetime
import argparse
import matplotlib.pyplot as plt

from pymatflow.qe.post.scf import scf_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of static scf running", type=str, default="tmp-qe-static")
    parser.add_argument("-f", "--file", help="output of static scf running: name only without path", type=str, default="static-scf.out")

    args = parser.parse_args()

    os.chdir(args.directory)
    os.system("mkdir -p post-processing")
    os.chdir("post-processing")

    scf = scf_out()
    scf.get_info(file="../%s" % args.file)

    # now in the directory of post-processing, output informations:

    # plot the scf energies
    plt.plot(scf.run_info["scf_energies"])
    plt.title("Total energy per scf step")
    plt.xlabel("Scf step")
    plt.ylabel("Total energy")
    plt.tight_layout()
    plt.savefig("total-energy-per-scf-step.png")
    plt.close()

    with open("scf-info.md", 'w', encoding='utf-8') as fout:
        fout.write("# SCF实验统计\n")
        fout.write("## SCF参数\n")
        for item in scf.scf_params:
            fout.write("- %s: %s\n" % (item, str(scf.scf_params[item])))
        fout.write("## 运行信息\n")
        # calculate the running time and print it out
        # depending on the value of seconds in the time string, there are three situations:
        # when second is smaller than 10  it will be divided from xx:xx x, and when second
        # is larger or equal to 10 it will be together(xx:xx:xx).
        # when all the three xx is smaller than 1 it will be like x: x: x
        # so we have to preprocess it to build the right time string to pass into
        # datetime.datetime.strptime()
        if len(scf.run_info["start_time"].split()) == 8:
            start_str = scf.run_info["start_time"].split()[5]+"-"+scf.run_info["start_time"].split()[7]
        elif len(scf.run_info["start_time"].split()) == 9:
            start_str = scf.run_info["start_time"].split()[5]+"-"+scf.run_info["start_time"].split()[7]+scf.run_info["start_time"].split()[8]
        elif len(scf.run_info["start_time"].split()) == 10:
            start_str = scf.run_info["start_time"].split()[5]+"-"+scf.run_info["start_time"].split()[7]+scf.run_info["start_time"].split()[8]+scf.run_info["start_time"].split()[9]
        else:
            print("===============================================\n")
            print("                  Warning !!!\n")
            print("===============================================\n")
            print("qe.post.scf.markdown_report:\n")
            print("failed to parse start_time string\n")
            sys.exit(1)
        if len(scf.run_info["stop_time"].split()) == 7:
            stop_str = scf.run_info["stop_time"].split()[6]+"-"+scf.run_info["stop_time"].split()[5]
        elif len(scf.run_info["stop_time"].split()) == 8:
            stop_str = scf.run_info["stop_time"].split()[7]+"-"+scf.run_info["stop_time"].split()[5]+scf.run_info["stop_time"].split()[6]
        elif len(scf.run_info["stop_time"].split()) == 9:
            stop_str = scf.run_info["stop_time"].split()[8]+"-"+scf.run_info["stop_time"].split()[5]+scf.run_info["stop_time"].split()[6]+scf.run_info["stop_time"].split()[7]
        else:
            print("===============================================\n")
            print("                  Warning !!!\n")
            print("===============================================\n")
            print("qe.post.scf.markdown_report:\n")
            print("failed to parse stop_time string\n")
            sys.exit(1)
        start = datetime.datetime.strptime(start_str, "%d%b%Y-%H:%M:%S")
        stop = datetime.datetime.strptime(stop_str, "%d%b%Y-%H:%M:%S")
        delta_t = stop -start
        fout.write("- Time consuming:\n")
        fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
        # end the time information
        for item in scf.run_info:
            fout.write("- %s: %s\n" % (item, str(scf.run_info[item])))

        fout.write("## 运行信息图示\n")
        fout.write("Total energy per scf step\n")
        fout.write("![Total energy per scf step](total-energy-per-scf-step.png)\n")

    # end of information output
    os.chdir("../")
    os.chdir("../")


    # --------------------------------------------------------------------------
    # print information to the terminal
    # --------------------------------------------------------------------------
    print("=====================================================================\n")
    print("                           post-qe-scf.py\n")
    print("---------------------------------------------------------------------\n")
    print("\n")
