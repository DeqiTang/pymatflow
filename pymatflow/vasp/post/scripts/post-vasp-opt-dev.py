#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import datetime
import argparse
from pymatflow.vasp.post.opt import opt_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of opt running", type=str, default="tmp-vasp-opt")

    args = parser.parse_args()

    os.chdir(args.directory)
    opt = opt_out()
    opt.get_info(outcar="OUTCAR", poscar="POSCAR")

    os.system("mkdir -p post-processing")
    os.chdir("post-processing")

    # now we are in post-processing, generate the output and return
    opt.print_trajectory()
    opt.print_final_structure()

    #plt.plot(opt.run_info["opt-energies"])
    #plt.title("Energy per scf step")
    #plt.xlabel("Scf step")
    #plt.ylabel("Total energy")
    #plt.tight_layout()
    #plt.savefig("energy-per-scf-step.png")
    #plt.close()

    with open("opt-info.md", 'w', encoding='utf-8') as fout:
        fout.write("# 几何优化实验统计\n")
        fout.write("几何优化类型: ISIF = %d\n" % opt.run_params["ISIF"])
        fout.write("几何优化任务是否结束:%s\n" % str(opt.job_done))
        if opt.job_done == True:
            fout.write("是否成功优化: %s\n" % str(opt.relaxed))
        else:
            fout.write("是否成功优化: %s\n" % ("运行未结束, 结果未知"))
        fout.write("## 离子步参数\n")
        for item in opt.run_params:
            fout.write("- %s: %s\n" % (item, str(opt.run_params[item])))
        fout.write("## 电子步参数\n")
        for item in opt.run_params:
            fout.write("- %s: %s\n" % (item, str(opt.run_params[item])))
        fout.write("## 运行信息\n")
        # calculate the running time and print it out
        # Importante: the length of the time string might be different, depending
        # on the value of hours and minutes and seconds. if they are two digits
        # number, they will be divided like: '11: 6: 2', only when they all are
        # two digtis number, they will not be divided '11:16:12'
        # so we have to preprocess it to build the right time string to pass into
        # datetime.datetime.strptime()
        start_str = opt.run_info["start_time"].split()[4]+"-"+opt.run_info["start_time"].split()[5]
        if opt.job_done == True:
            #stop_str = opt.run_info["stop-time"].split()[8]+"-"+opt.run_info["stop-time"].split()[5]+opt.run_info["stop-time"].split()[6]+opt.run_info["stop-time"].split()[7]
            pass

        start = datetime.datetime.strptime(start_str, "%Y.%m.%d-%H:%M:%S")
        #if opt.job_done == True:
        #    stop = datetime.datetime.strptime(stop_str, "%d%b%Y-%H:%M:%S")
        #    delta_t = stop -start
        fout.write("- Time consuming:\n")
        fout.write("  - job starts at %s\n" % start)
        fout.write("  - Elapsed time: %.3f(sec) = %.3f(min) = %.3f(hour)\n" % (opt.run_info["elapsed_time"], opt.run_info["elapsed_time"]/60, opt.run_info["elapsed_time"]/3600))
        #if opt.job_done == True:
        #    fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
        #else:
        #    fout.write("  - job is not finished yet, but it starts at %s\n" % start)
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
        fout.write("![Total forces per SCF](total-forces-rms-per-scf.png)\n")


    os.chdir("../")
    os.chdir("../")


    # --------------------------------------------------------------------------
    # print information to the terminal
    # --------------------------------------------------------------------------
    print("=====================================================================\n")
    print("                           post-vasp-opt-dev.py\n")
    print("---------------------------------------------------------------------\n")
    print("\n")
