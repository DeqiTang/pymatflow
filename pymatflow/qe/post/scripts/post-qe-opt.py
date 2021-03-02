#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import datetime
import argparse
import matplotlib.pyplot as plt

from pymatflow.qe.post.opt import OptOut


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", type=str, help="directory of optimization running")
    parser.add_argument("-f", "--file", type=str, help="output of optimization running: name only without path")

    args = parser.parse_args()

    os.chdir(args.directory)
    os.system("mkdir -p post-processing")
    os.chdir("post-processing")

    opt = OptOut()
    opt.get_info(file="../%s" % args.file)

    # now in the directory of post-processing, output informations:

    #plt.plot(opt.run_info["iterations"])
    #plt.title("Iterations per SCF")
    #plt.xlabel("Scf cycles")
    #plt.ylabel("iterations")
    #plt.tight_layout()
    #plt.savefig("iterations-per-scf.png")
    #plt.close()

    #plt.plot(opt.run_info["total-energies"])
    #plt.title("Total energies per SCF")
    #plt.xlabel("Scf cycles")
    #plt.ylabel("Total Energies (Ry)")
    #plt.tight_layout()
    #plt.savefig("total-energies-per-scf.png")
    #plt.close()

    #plt.plot(opt.run_info["fermi-energies"])
    #plt.title("Fermi energies per SCF")
    #plt.xlabel("Scf cycles")
    #plt.ylabel("Fermi energies (eV)")
    #plt.tight_layout()
    #plt.savefig("fermi-energies-per-scf.png")
    #plt.close()

    #plt.plot(opt.run_info["total-forces"])
    #plt.title("Total forces per SCF")
    #plt.xlabel("Scf cycles")
    #plt.ylabel("Total forces (Ry/au)")
    #plt.tight_layout()
    #plt.savefig("total-forces-per-scf.png")
    #plt.close()

    with open("optimization-info.md", 'w', encoding='utf-8') as fout:
        fout.write("# 几何优化实验统计\n")
        fout.write("几何优化类型: %s\n" % opt.run_type)
        fout.write("几何优化任务是否结束:%s\n" % str(opt.job_done))
        if opt.job_done == True:
            fout.write("是否成功优化: %s\n" % str(opt.relaxed))
        else:
            fout.write("是否成功优化: %s\n" % ("运行未结束, 结果未知"))
        fout.write("## 优化参数\n")
        for item in opt.opt_params:
            fout.write("- %s: %s\n" % (item, str(opt.opt_params[item])))
        fout.write("## 运行信息\n")
        # calculate the running time and print it out
        # Importante: the length of the time string might be different, depending
        # on the value of hours and minutes and seconds. if they are two digits
        # number, they will be divided like: '11: 6: 2', only when they all are
        # two digtis number, they will not be divided '11:16:12'
        # so we have to preprocess it to build the right time string to pass into
        # datetime.datetime.strptime()
        if len(opt.run_info["start_time"].split()) == 8:
            start_str = opt.run_info["start_time"].split()[5]+"-"+opt.run_info["start_time"].split()[7]
        elif len(opt.run_info["start_time"].split()) == 9:
            start_str = opt.run_info["start_time"].split()[5]+"-"+opt.run_info["start_time"].split()[7]+opt.run_info["start_time"].split()[8]
        elif len(opt.run_info["start_time"].split()) == 10:
            start_str = opt.run_info["start_time"].split()[5]+"-"+opt.run_info["start_time"].split()[7]+opt.run_info["start_time"].split()[8]+opt.run_info["start_time"].split()[9]
        else:
            print("===============================================\n")
            print("                  Warning !!!\n")
            print("===============================================\n")
            print("post-qe-opt.py:\n")
            print("failed to parse start_time string\n")
            sys.exit(1)
        if opt.job_done == True:
            if len(opt.run_info["stop_time"].split()) == 7:
                stop_str = opt.run_info["stop_time"].split()[6]+"-"+opt.run_info["stop_time"].split()[5]
            elif len(opt.run_info["stop_time"].split()) == 8:
                stop_str = opt.run_info["stop_time"].split()[7]+"-"+opt.run_info["stop_time"].split()[5]+opt.run_info["stop_time"].split()[6]
            elif len(opt.run_info["stop_time"].split()) == 9:
                stop_str = opt.run_info["stop_time"].split()[8]+"-"+opt.run_info["stop_time"].split()[5]+opt.run_info["stop_time"].split()[6]+opt.run_info["stop_time"].split()[7]
            else:
                print("===============================================\n")
                print("                  Warning !!!\n")
                print("===============================================\n")
                print("post-qe-opt.pyt:\n")
                print("failed to parse stop_time string\n")
                sys.exit(1)

        start = datetime.datetime.strptime(start_str, "%d%b%Y-%H:%M:%S")
        if opt.job_done == True:
            stop = datetime.datetime.strptime(stop_str, "%d%b%Y-%H:%M:%S")
            delta_t = stop -start
        fout.write("- Time consuming:\n")
        if opt.job_done == True:
            fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
        else:
            fout.write("  - job is not finished yet, but it starts at %s\n" % start)
        # end the time information
        for item in opt.run_info:
            fout.write("- %s: %s\n" % (item, str(opt.run_info[item])))

        fout.write("## 运行信息图示\n")
        fout.write("Iterations per SCF\n")
        fout.write("![Iterations per SCF](iterations-per-scf.png)\n")

        fout.write("Total energies per SCF\n")
        fout.write("![Total energies per SCF](total-energies-per-scf.png)\n")

        fout.write("Fermi energies per SCF\n")
        fout.write("![Fermi energies per SCF](fermi-energies-per-scf.png)\n")

        fout.write("Total forces per SCF\n")
        fout.write("![Total forces per SCF](total-forces-per-scf.png)\n")


    os.chdir("../")
    os.chdir("../")


    # --------------------------------------------------------------------------
    # print information to the terminal
    # --------------------------------------------------------------------------
    print("=====================================================================\n")
    print("                           post-qe-opt.py\n")
    print("---------------------------------------------------------------------\n")
    print("\n")
