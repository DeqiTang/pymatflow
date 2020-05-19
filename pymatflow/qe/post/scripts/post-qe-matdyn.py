#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse


"""
usage:
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously projwfc running directory", type=str, default="tmp-qe-static")

    parser.add_argument("--matdyn-in", help="input file for the matdyn.x calculation", type=str, default="matdyn.in")

    parser.add_argument("--option", type=str, default="gnuplot",
            choices=["matplotlib", "gnuplot"],
            help="choosing gnuplot or matplotlib to do the band plot")

    parser.add_argument("--freq", type=float, nargs="+",
            default=[0, 1.0],
            help="band range to plot. in percentage")

    args = parser.parse_args()


    # get the x coordinate and label for high symmetry k point from matdyn.x input file
    with open(os.path.join(args.directory, args.matdyn_in), 'r') as fin:
        matdyn_in = fin.readlines()

    specialk = []
    # specialk: [{'label': 'GAMMA', 'xcoord': float}, ...]
    for line in matdyn_in:
        if len(line.split("#")) == 0:
            continue
        if len(line.split()) > 4 and len(line.split("#")[0].split()) == 4:
            print(line)
            specialk.append({
                "label": line.split("\n")[0].split("#")[1].upper(),
                "xcoord": float(line.split()[3]),
                })
    

    locs = [specialk[i]["xcoord"] for i in range(len(specialk))]
    labels = ["{/symbol G}" if specialk[i]["label"] == "GAMMA" else "%s" % specialk[i]["label"] for i in range(len(specialk))]
    #
    # sometime the xcoord of two specialk might be the same
    # either caused from physical reason or when you specif
    # 0 to connect the two special k point.
    # whatever the reason we should join the two label to one
    # label like K|U
    locs_refined = []
    labels_refined = []
    labels_refined.append(labels[0])
    locs_refined.append(locs[0])
    for i in range(1, len(labels)):
        if locs[i]  == locs[i-1]:
            # join labels[i] and labels[i-1]
            labels_refined[-1] = "%s | %s" % (labels[i-1], labels[i])
        else:
            labels_refined.append(labels[i])
            locs_refined.append(locs[i])


    # get flfrq file name
    for line in matdyn_in:
        if len(line.split("=")) == 0:
            continue
        if line.split("=")[0].split()[0] == 'flfrq':
            flfrq = line.split("=")[1].split("\n")[0].split()[0]
            break


    os.chdir(args.directory)

    # generating bash script file to get the gnuplot 
    with open("matdyn-analysis.sh", 'w') as fout:
        fout.write("#!/bin/bash\n")
        fout.write("\n")
        fout.write("flfrq=%s\n" % flfrq)
        fout.write("nbnd=`cat ${flfrq} | head -n +1 | cut -d \",\" -f 1 | cut -d \" \" -f 5`\n")
        fout.write("range_min=%.2f\n" % args.freq[0])
        fout.write("range_max=%.2f\n" % args.freq[1])
        fout.write("column_begin=`echo 2 $nbnd ${range_min} | awk '{ printf \"%d\", $1+$2*$3 }'`\n")
        fout.write("column_end=`echo 2 $nbnd ${range_max} | awk '{ printf \"%d\", $1+$2*$3-1 }'`\n")
        fout.write("cat > plot.gnuplot<<EOF\n")
        fout.write("set terminal gif\n")
        fout.write("set output 'phonon-band.gif'\n")
        fout.write("unset key\n")
        fout.write("set parametric\n")
        fout.write("set title 'Phonon spectrum'\n")
        fout.write("set xlabel 'K'\n")
        fout.write("set ylabel 'Frequency (cm^{-1})'\n")

        fout.write("set xtics(")
        for i in range(int(len(labels_refined)-1)):
            fout.write("'%s' %f, " % (labels_refined[i], locs_refined[i]))
        fout.write("'%s' %f)\n" % (labels_refined[-1], locs_refined[-1]))

        fout.write("set grid xtics ytics\n")
        fout.write("set autoscale\n")
        fout.write("plot for [i=${column_begin}:${column_end}:1] 'matdyn.freq.gp' using 1:i w l\n")
        fout.write("EOF\n")
        fout.write("gnuplot plot.gnuplot\n")

    os.system("bash matdyn-analysis.sh")
    os.chdir("../")
