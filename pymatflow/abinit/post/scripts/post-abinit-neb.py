#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import argparse

#from pymatflow.abinit.post.neb import neb_post
from pymatflow.abinit.post.neb import neb_out



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously neb running directory", type=str, default="tmp-abinit-neb")
    parser.add_argument("--nebout", help="output file of neb calculation", type=str, default="neb.out")

    args = parser.parse_args()

    #os.chdir(args.directory)
    #task = neb_post(args.nebout)
    #task.export()
    #os.chdir("../")


    os.chdir(args.directory)
    neb = neb_out()
    neb.get_info(file=args.enbout)
    os.system("mkdir -p post-processing")

    #print_trajectory
    with open("trajectory_initial.xyz", 'w') as fout:
        for i in range(len(neb.trajectory_initial)):
            fout.write("%d\n" % len(neb.trajectory_initial[i]))
            fout.write("i = %d\n" % i)
            for atom in neb.trajectory_initial[i]:
                fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
    with open("trajectory_final.xyz", 'w') as fout:
        for i in range(len(neb.trajectory_final)):
            fout.write("%d\n" % len(neb.trajectory_final[i]))
            fout.write("i = %d\n" % i)
            for atom in neb.trajectory_final[i]:
                fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))

#def view_trajectory(neb, trajfile_initial="trajectory-initial.xyz", trajfile_final="trajectory-final.xyz"):
    #os.system("xcrysden --xyz %s" % trajfile)
    #subprocess.call(["xcrysden", "--xyz", trajfile_initial])
    #subprocess.call(["xcrysden", "--xyz", trajfile_final])

    #plot_run_info
    plt.plot([neb.run_info["etotal-per-image"][i][1] for i in range(len(neb.run_info["etotal-per-image"]))])
    plt.title("Total energies per image")
    plt.xlabel("image")
    plt.ylabel("Total Energies (Hartree)")
    plt.tight_layout()
    plt.savefig("etotal-per-image.png")
    plt.close()



    with open("neb-info.md", 'w', encoding='utf-8') as fout:
        fout.write("# 过渡态搜索优化实验统计\n")
        fout.write("## 过渡态参数\n")
        for item in neb.neb_params:
            fout.write("- %s: %s\n" % (item, str(neb.neb_params[item])))
        fout.write("## 运行信息\n")
        # calculate the running time and print it out

        # end the time information
        for item in neb.run_info:
            fout.write("- %s: %s\n" % (item, str(neb.run_info[item])))

        fout.write("## 运行信息图示\n")

        fout.write("Total energies per image\n")
        fout.write("![Total energies per image](etotal-per-image.png)\n")

    os.chdir("../")
    os.chdir("../")


    # --------------------------------------------------------------------------
    # print information to the terminal
    # --------------------------------------------------------------------------
    print("=====================================================================\n")
    print("                           post-abinit-neb.py\n")
    print("---------------------------------------------------------------------\n")
    print("\n")
