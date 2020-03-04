#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

import numpy as np
import matplotlib.pyplot as plt


"""
该脚本在左图中定位高对称点的位置的方法是，读取`--kpath`文件中定义的
高对称点信息，根据其长度以及"|"的数量(n_div)可以知道bandline的
数量(n_band_line):`len(kpath) - 1 -n_div`。知道了bandline的数量我们
就可以从能带横坐标数据的长度(`len(xcoord)`)来获取bandline的长度(有
多少个x点):`len(xcoord) / n_band_line`这个值实际上等于KPOINTS中设置
的链接高对称点的数值。然后根据这些信息，我们就可以把高对称点的符号
与xcoord联系在一起。
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--kpath", type=str, nargs="+", default=None,
            help="manual input kpath for band structure calculation")

    parser.add_argument("--kpath-file", type=str, default="kpath-from-seekpath.txt",
            help="file to read the kpath for band structure calculation")

    parser.add_argument("--bands-p4vasp", type=str, default=None,
            help="bands data exported by p4vasp in gnuplot format")

    parser.add_argument("--bandrange", type=float, nargs="+",
            default=[0, 1.0],
            help="band range to plot. in percentage")

    parser.add_argument("--engine", type=str, default="matplotlib",
            choices=["matplotlib", "gnuplot"],
            help="plot engine, matplotlib or gnuplot")

    parser.add_argument("--output-dir", type=str, default="./",
           help="directory to put the analysis result files")


    args = parser.parse_args()


    if args.kpath != None:
        # kpath from script argument args.kpath
        kpath = []
        for kpoint in args.kpath:
            if kpoint.split()[4] != "|":
                kpath.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    int(kpoint.split()[4]),
                    ])
            elif kpoint.split()[4] == "|":
                kpath.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    "|",
                    ])
    elif args.kpath == None:
        # kpath read from file specified by args.kpath_file
        # file is in format like this
        """
        5
        0.0 0.0 0.0 #GAMMA 15
        x.x x.x x.x #XXX |
        x.x x.x x.x #XXX 10
        x.x x.x x.x #XXX 15
        x.x x.x x.x #XXX 20
        """
        # if there is a '|' behind the label it means the path is
        # broken after that point!!!
        kpath = []
        with open(args.kpath_file, 'r') as fin:
            kpath_file = fin.readlines()
        nk = int(kpath_file[0])
        for i in range(nk):
            if kpath_file[i+1].split("\n")[0].split()[4] != "|":
                kpath.append([
                    float(kpath_file[i+1].split()[0]),
                    float(kpath_file[i+1].split()[1]),
                    float(kpath_file[i+1].split()[2]),
                    kpath_file[i+1].split()[3].split("#")[1].upper(),
                    int(kpath_file[i+1].split()[4]),
                    ])
            elif kpath_file[i+1].split("\n")[0].split()[4] == "|":
                kpath.append([
                    float(kpath_file[i+1].split()[0]),
                    float(kpath_file[i+1].split()[1]),
                    float(kpath_file[i+1].split()[2]),
                    kpath_file[i+1].split()[3].split("#")[1].upper(),
                    '|',
                    ])
    else:
        pass
        # 2 not in args.printout_option
        # do not calculate the band structure
        # no need to set the kpath


    #
    with open(args.bands_p4vasp, 'r') as fin:
        bands_data = fin.readlines()

    xcoord = []
    for line in bands_data:
        if line == '\n':
            break
        xcoord.append(float(line.split()[0]))

    bandsenergy = [] # [np.array([x, x, x, x, x, x, x, x....]), np.array([x, x, x, x, x, x, x....]), .......]
    dataall = np.loadtxt(bands_data)
    for i in range(int(len(dataall) / len(xcoord))):
        bandsenergy.append(dataall[i*len(xcoord) : (i+1)*len(xcoord), 1])


    locs = []
    labels_for_matplotlib = []
    labels_for_gnuplot = []
    n_div = 0
    for kpoint in kpath:
        if kpoint[4] == "|":
            n_div += 1
    n_band_line = len(kpath) - 1 - n_div
    intersections = len(xcoord) / n_band_line
    for i in range(n_band_line):
        print(intersections*i)
        locs.append(xcoord[int(intersections*i)])
    locs.append(xcoord[-1])

    labels_for_matplotlib.append(r"$%s$" % kpath[0][3].upper() if kpath[0][3].upper() != "GAMMA" else r"$\Gamma$")
    labels_for_gnuplot.append("%s" % kpath[0][3].upper() if kpath[0][3].upper() != "GAMMA" else "{/symbol G}")
    for i in range(1, len(kpath)):
        if kpath[i-1][4] != "|":
            labels_for_matplotlib.append(r"$%s$" % kpath[i][3].upper() if kpath[i][3].upper() != "GAMMA" else r"$\Gamma$")
            labels_for_gnuplot.append("%s" % kpath[i][3].upper() if kpath[i][3].upper() != "GAMMA" else "{/symbol G}")
        else:
            labels_for_matplotlib[-1] = r"$%s | %s$" % (labels_for_matplotlib[-1].split("$")[1], kpath[i][3].upper())
            labels_for_gnuplot[-1] = "%s | %s" % (labels_for_gnuplot[-1], kpath[i][3].upper())


    begin = int(len(bandsenergy) * args.bandrange[0])
    end = int(len(bandsenergy) * args.bandrange[1])

    if args.engine == "matplotlib":
        for i in range(begin, end):
            plt.plot(xcoord, bandsenergy[i])
        plt.title("Band Structure")
        plt.ylabel(r"$\mathit{E}-\mathit{E}_\mathrm{f} (eV)$")
        plt.xlabel("Kpoints")
        plt.grid(which="major", axis="x", linewidth=0.75, linestyle="-", color="0.75")
        plt.grid(which="major", axis="y", linewidth=0.75, linestyle="-", color="0.75")
        plt.tight_layout()
        plt.xticks(locs, labels_for_matplotlib)
        plt.savefig(os.path.join(args.output_dir, "bandstructure.png"))
        plt.show()

    if args.engine == "gnuplot":
        # output the band structure data to band-structure.data.gnu which contains part or all
        # of the energy data depending on args.bandrange
        with open(os.path.join(args.output_dir, "band-structure.data.gnu"), 'w') as fout:
            for i in range(begin, end):
                for j in range(len(xcoord)):
                    fout.write("%f %f\n" % (xcoord[j], bandsenergy[i][j]))
                fout.write("\n")
        with open(os.path.join(args.output_dir, "band-structure-plot.gnuplot"), 'w') as fout:
            fout.write("set terminal gif\n")
            fout.write("set output '%s'\n" % os.path.join(args.output_dir, "band-structure.gif"))
            fout.write("unset key\n")
            fout.write("set parametric\n")
            fout.write("set title 'Band Structure'\n")
            fout.write("set xlabel 'K'\n")
            fout.write("set ylabel 'Energy(eV)'\n")
            fout.write("set xtics(")
            for i in range(len(labels_for_gnuplot)):
                fout.write("'%s' %f, " % (labels_for_gnuplot[i], locs[i]))
            fout.write(")\n")
            fout.write("set grid xtics ytics\n")
            fout.write("set autoscale\n")
            fout.write("plot '%s' w l\n" % os.path.join(args.output_dir, "band-structure.data.gnu"))
        os.system("cd %s; gnuplot band-structure-plot.gnuplot" % args.output_dir)
        os.system("eog %s\n" % os.path.join(args.output_dir, "band-structure.gif"))
