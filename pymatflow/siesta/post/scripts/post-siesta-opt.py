#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
#from pymatflow.siesta.post.opt import opt_post
from pymatflow.siesta.post.opt import opt_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of geometric optimization running", type=str, default="tmp-siesta-opt")
    parser.add_argument("-f", "--file", help="output of opt running", type=str, default="geometric-optimization.out")
    parser.add_argument("--view-traj", type=str, default="yes",
            choices=["yes", "no"],
            help="whether to view the trajectory.")

    args = parser.parse_args()

    #os.chdir(args.directory)
    #task = opt_post(outputfile=args.file)
    #task.export()
    #if args.view_traj == "yes":
    #    task.view_trajectory()
    #os.chdir("../")


    os.chdir("args.directory")
    opt = opt_out()
    opt.get_info(file=args.file)
    os.system("mkdir -p post-processing")

    #view_trajectory # trajfile="siesta.ANI"):
    #trajfile = "siesta.ANI"
    #subprocess.call(["xcrysden", "--xyz", trajfile])

    #plot_run_info
    plt.plot(opt.run_info["total-energies"])
    plt.title("Total energies per SCF")
    plt.xlabel("Scf cycles")
    plt.ylabel("Total Energies (eV)")
    plt.tight_layout()
    plt.savefig("total-energies-per-scf.png")
    plt.close()

    plt.plot(opt.run_info["scf-iterations-converged"])
    plt.title("Iterations per SCF(converged)")
    plt.xlabel("Scf cycle(converged)")
    plt.ylabel("Scf iterations")
    plt.tight_layout()
    plt.savefig("iterations-per-scf-converged.png")
    plt.close()


    """
    when writing Chinese to a file you must specify
    encoding='utf-8' when open the file for writing
    """
    with open("opt-info.md", 'w', encoding='utf-8') as fout:
        fout.write("# 几何优化实验统计\n")
        fout.write("是否进行变胞优化: %s\n" % opt.opt_params["VariableCell"])
        fout.write("几何优化任务是否结束:%s\n" % str(opt.job_completed))
        if opt.job_completed == True:
            fout.write("是否成功优化: %s\n" % str(opt.relaxed))
        else:
            fout.write("是否成功优化: %s\n" % "运行未结束, 结果未知")
        fout.write("## 优化参数\n")
        for item in opt.opt_params:
            fout.write("- %s: %s\n" % (item, str(opt.opt_params[item])))
        fout.write("## 运行信息\n")
        # calculate the running time and print it out
        start = datetime.datetime.strptime(opt.run_info["start-time"].split()[4]+"-"+opt.run_info["start-time"].split()[5], "%d-%b-%Y-%H:%M:%S")
        if opt.job_completed == True:
            stop = datetime.datetime.strptime(opt.run_info["stop-time"].split()[4]+"-"+opt.run_info["stop-time"].split()[5], "%d-%b-%Y-%H:%M:%S")
            delta_t = stop -start
        fout.write("- Time consuming:\n")
        if opt.job_completed == True:
            fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
        else:
            fout.write("  - job is not finished yet, but it starts at %s\n" % start)
        # end the time information
        for item in opt.run_info:
            fout.write("- %s: %s\n" % (item, str(opt.run_info[item])))

        fout.write("## 运行信息图示\n")
        fout.write("Iterations per SCF\n")
        fout.write("![Iterations per SCF](iterations-per-scf-converged.png)\n")

        fout.write("Total energies per SCF\n")
        fout.write("![Total energies per SCF](total-energies-per-scf.png)\n")

    os.chdir("../")
    os.chdir("../")
