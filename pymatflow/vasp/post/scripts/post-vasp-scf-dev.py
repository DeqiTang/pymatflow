#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import datetime
import argparse
from pymatflow.vasp.post.scf import scf_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of scf running", type=str, default="tmp-vasp-static")

    args = parser.parse_args()

    os.chdir(args.directory)
    scf = scf_out()
    scf.get_info(outcar="OUTCAR")

    os.system("mkdir -p post-processing")
    os.chdir("post-processing")

    # now we are in post-processing, generate the output and return

    plt.plot(self.run_info["scf-energies"])
    plt.title("Energy per scf step")
    plt.xlabel("Scf step")
    plt.ylabel("Total energy")
    plt.tight_layout()
    plt.savefig("energy-per-scf-step.png")
    plt.close()

    with open("scf-info.md", 'w', encoding='utf-8') as fout:
        fout.write("# SCF实验统计\n")
        fout.write("## SCF参数\n")
        for item in self.scf_params:
            fout.write("- %s: %s\n" % (item, str(self.scf_params[item])))
        fout.write("## 运行信息\n")
        # calculate the running time and print it out
        # Importante: the length of the time string might be different, depending
        # on the value of hours and minutes and seconds. if they are two digits
        # number, they will be divided like: '11: 6: 2', only when they all are
        # two digtis number, they will not be divided '11:16:12'
        # so we have to preprocess it to build the right time string to pass into
        # datetime.datetime.strptime()
        start_str = self.run_info["start-time"].split()[4]+"-"+self.run_info["start-time"].split()[5]
        if self.job_done == True:
            #stop_str = self.run_info["stop-time"].split()[8]+"-"+self.run_info["stop-time"].split()[5]+self.run_info["stop-time"].split()[6]+self.run_info["stop-time"].split()[7]
            pass

        start = datetime.datetime.strptime(start_str, "%Y.%m.%d-%H:%M:%S")
        #if self.job_done == True:
        #    stop = datetime.datetime.strptime(stop_str, "%d%b%Y-%H:%M:%S")
        #    delta_t = stop -start
        fout.write("- Time consuming:\n")
        fout.write("  - job starts at %s\n" % start)
        fout.write("  - Elapsed time: %.3f(sec) = %.3f(min) = %.3f(hour)\n" % (self.run_info["elapsed-time"], self.run_info["elapsed-time"]/60, self.run_info["elapsed-time"]/3600))
        #if self.job_done == True:
        #    fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
        #else:
        #    fout.write("  - job is not finished yet, but it starts at %s\n" % start)
        # end the time information
        for item in self.run_info:
            fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

        fout.write("## 运行信息图示\n")
        fout.write("Total energy per scf step\n")
        fout.write("![energy per scf step](energy-per-scf-step.png)\n")


    os.chdir("../")
    os.chdir("../")


    # --------------------------------------------------------------------------
    # print information to the terminal
    # --------------------------------------------------------------------------
    print("=====================================================================\n")
    print("                           post-vasp-scf-dev.py\n")
    print("---------------------------------------------------------------------\n")
    print("\n")
