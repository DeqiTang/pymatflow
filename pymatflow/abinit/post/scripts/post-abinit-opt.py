#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.abinit.post.opt import opt_out



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously optimization running directory", type=str, default="tmp-abinit-opt")
    parser.add_argument("--optout", help="output file of optimization calculation", type=str, default="geometric-optimization.out")

    args = parser.parse_args()

    #os.chdir(args.directory)
    #task = opt_post(args.optout)
    #task.export()
    #os.chdir("../")

    os.chdir(args.directory)
    opt = opt_out()
    opt.get_info(file=args.output)

    os.system("mkdir -p post-processing")
    #


    #print_trajectory
    with open("trajectory.xyz", 'w') as fout:
        for i in range(len(opt.trajectory)):
            fout.write("%d\n" % len(opt.trajectory[i]))
            fout.write("i = %d\n" % i)
            for atom in opt.trajectory[i]:
                fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))

    #subprocess.call(["xcrysden", "--xyz", trajfile])

    #print_final_structure
    bohr = 0.5291772108
    if len(opt.acells) == 0 and len(opt.rprimds) == 0:
        # not optimized cell
        cell = copy.deepcopy(opt.outvars_before["rprim"])
        for n in range(3):
            cell[n] = cell[n] * opt.outvars_before["acell"][0] * bohr
            cell[n+3] = cell[n+3] * opt.outvars_before["acell"][1] * bohr
            cell[n+6] = cell[n+6] * opt.outvars_before["acell"][2] * bohr
        #
    else:
        cell = opt.cells[-1]
    #
    with open("final-structure.xyz", 'w') as fout:
        fout.write("%d\n" % len(opt.trajectory[-1]))
        fout.write("cell: %.9f %.9f %.9f | %.9f %.9f %.9f | %.9f %.9f %.9f\n" % (cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], cell[6], cell[7], cell[8]))
        for atom in opt.trajectory[-1]:
            fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))

    #plot_run_info
    plt.plot(opt.run_info["iterations"])
    plt.title("Iterations per SCF")
    plt.xlabel("Scf cycles ")
    plt.ylabel("iterations")
    plt.tight_layout()
    plt.savefig("iterations-per-scf.png")
    plt.close()

    plt.plot(opt.run_info["total-energies"])
    plt.title("Total energies per SCF")
    plt.xlabel("Scf cycles")
    plt.ylabel("Total Energies (Hartree)")
    plt.tight_layout()
    plt.savefig("total-energies-per-scf.png")
    plt.close()


    with open("opt-info.md", 'w', encoding='utf-8') as fout:
        fout.write("# 几何优化实验统计\n")
        fout.write("## 几何优化参数\n")
        for item in opt.opt_params:
            fout.write("- %s: %s\n" % (item, str(opt.opt_params[item])))
        fout.write("## 运行信息\n")
        # calculate the running time and print it out

        # end the time information
        for item in opt.run_info:
            fout.write("- %s: %s\n" % (item, str(opt.run_info[item])))

        fout.write("## 运行信息图示\n")

        fout.write("Iterations per SCF\n")
        fout.write("![Iterations per SCF](iterations-per-scf.png)\n")

        fout.write("Total energies per scf\n")
        fout.write("![Total energies per SCF](total-energies-per-scf.png)\n")


    os.chdir("../")
    os.chdir('../')

    # --------------------------------------------------------------------------
    # print information to the terminal
    # --------------------------------------------------------------------------
    print("=====================================================================\n")
    print("                           post-abinit-opt.py\n")
    print("---------------------------------------------------------------------\n")
    print("\n")
