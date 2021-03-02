#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.base.xyz import BaseXyz

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of phonopy running", type=str, default="tmp-qe-phonopy")
    parser.add_argument("-f", "--file", help="input structure file", type=str, default=None)


    parser.add_argument("--qpath", type=str, nargs="+", default=None,
            help="manual input qpath with labels to set the BAND for phonopy analysis")

    parser.add_argument("--qpath-file", type=str, default="kpath-from-seekpath.txt",
            help="file to read the qpath to set the BAND for phonopy analysis")

    parser.add_argument("--mp", type=int, nargs="+",
            default=[8, 8, 8],
            help="MP for mesh.conf and pdos.conf")

    parser.add_argument("--supercell-n", type=int, nargs="+",
            default=None,
            help="supercell_n used in phonopy calculation")

    parser.add_argument("--engine", type=str, default="matplotlib",
            choices=["matplotlib", "gnuplot"],
            help="plot engine, matplotlib or gnuplot")
            
    args = parser.parse_args()

    engine = args.engine
    
    xyz = BaseXyz()
    xyz.get_xyz(args.file)

    # obtain the qpath
    qpath = [] # [[kx, ky, kz, label, end_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', None], ...]
    # [[kx, ky, kz, label, connect_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', 15], ...]
    # if connect_indicator in a kpoint is an integer, then it will connect to the following point
    # through the number of kpoints defined by connect_indicator.
    # if connect_indicator in a kpoint is '|', then it will not connect to the following point,
    if args.qpath != None:
        # qpath from script argument args.qpath
        for qpoint in args.qpath:
            if qpoint.split()[4] != "|":
                qpath.append([
                    float(qpoint.split()[0]),
                    float(qpoint.split()[1]),
                    float(qpoint.split()[2]),
                    qpoint.split()[3].upper(),
                    None, # this value is actually not used by phonopy so we just set it to None if it is not "|"
                    ])
            elif qpoint.split()[4] == "|":
                qpath.append([
                    float(qpoint.split()[0]),
                    float(qpoint.split()[1]),
                    float(qpoint.split()[2]),
                    qpoint.split()[3].upper(),
                    "|",
                    ])
    else:
        # qpath read from file specified by args.qpath_file
        # file is in format like this
        """
        4
        0.0 0.0 0.0 #GAMMA 15
        x.x x.x x.x #XXX |
        x.x x.x x.x #XXX 15
        x.x x.x x.x #XXX 10
        """
        with open(args.qpath_file, 'r') as fin:
            qpath_file = fin.readlines()
        nq = int(qpath_file[0])
        for i in range(nq):
            if qpath_file[i+1].split("\n")[0].split()[4] != "|":
                qpath.append([
                    float(qpath_file[i+1].split()[0]),
                    float(qpath_file[i+1].split()[1]),
                    float(qpath_file[i+1].split()[2]),
                    qpath_file[i+1].split()[3].split("#")[1].upper(),
                    None, # this value is actually not used by phonopy so we just set it to None if it is not "|"
                    ])
            elif qpath_file[i+1].split("\n")[0].split()[4] == "|":
                qpath.append([
                    float(qpath_file[i+1].split()[0]),
                    float(qpath_file[i+1].split()[1]),
                    float(qpath_file[i+1].split()[2]),
                    qpath_file[i+1].split()[3].split("#")[1].upper(),
                    '|',
                    ])
    #

    # get the disps information
    os.chdir(args.directory)
    os.system("ls | grep 'supercell-' > pos.data")
    disps = []
    with open("pos.data", 'r') as fin:
        for line in fin:
            disps.append(line.split(".")[0].split("-")[1])
    os.chdir("../")

    os.chdir(args.directory)

    # generate the result analysis bash script and necessary config files

    with open("mesh.conf", 'w') as fout:
        fout.write("ATOM_NAME =")
        for element in xyz.specie_labels:
            fout.write(" %s" % element)
        fout.write("\n")
        fout.write("DIM = %d %d %d\n" % (args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))
        fout.write("MP = %d %d %d\n" % (args.mp[0], args.mp[1], args.mp[2]))



    with open("pdos.conf", 'w') as fout:
        fout.write("ATOM_NAME =")
        for element in xyz.specie_labels:
            fout.write(" %s" % element)
        fout.write("\n")
        fout.write("DIM = %d %d %d\n" % (args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))
        fout.write("MP = %d %d %d\n" % (args.mp[0], args.mp[1], args.mp[2]))
        fout.write("PDOS = 1 2, 3 4 5 5\n")

    with open("band.conf", 'w') as fout:
        fout.write("ATOM_NAME =")
        for element in xyz.specie_labels:
            fout.write(" %s" % element)
        fout.write("\n")
        # tke use of PRIMITIVE_AXES will find the primitive cell of the structure
        # and use it to analyse the phonon band structure
        # however, the use of primitive cell will not affect the q path setting
        # so whether we use PRIMITIVE cell or not, we can set the same q path
        fout.write("PRIMITIVE_AXES = AUTO\n") # we can also specify a matrix, but AUTO is recommended now in phonopy
        fout.write("GAMMA_CENTER = .TRUE.\n")
        fout.write("BAND_POINTS = 101\n")
        fout.write("BAND_CONNECTION = .TRUE.\n")
        fout.write("DIM = %d %d %d\n" % (args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))
        fout.write("BAND =")
        for i in range(len(qpath) - 1):
            if qpath[i][4] == None:
                fout.write(" %f %f %f" % (qpath[i][0], qpath[i][1], qpath[i][2]))
            elif qpath[i][4] == "|":
                fout.write(" %f %f %f," % (qpath[i][0], qpath[i][1], qpath[i][2]))
            else:
                pass
        fout.write(" %f %f %f" % (qpath[-1][0], qpath[-1][1], qpath[-1][2]))
        fout.write("\n")
        
        fout.write("BAND_LABELS =")
        for i in range(len(qpath) - 1):
            if qpath[i][4] == None:
                if qpath[i][3].upper() == "GAMMA":
                    fout.write(" $\Gamma$")
                else:
                    fout.write(" $%s$" % qpath[i][3])
            elif qpath[i][4] == "|":
                if qpath[i][3].upper() == "GAMMA":
                    fout.write(" $\Gamma$,")
                else:
                    fout.write(" $%s$," % qpath[i][3])
            else:
                pass
        if qpath[-1][3].upper() == "GAMMA":
            fout.write(" $\Gamma$")
        else:
            fout.write(" $%s$" % qpath[-1][3])            
        fout.write("\n")

    with open("phonopy-analysis.sh", 'w') as fout:
        fout.write("#!/bin/bash\n\n")
        fout.write("# generate the FORCE_SET\n")
        fout.write("phonopy --qe -f supercell-{001..%s}.out\n" % (disps[-1]))
        fout.write("# plot The density of states (DOS)\n")
        fout.write("phonopy --qe -p mesh.conf -c pos.in\n")
        fout.write("# Thermal properties are calculated with the sampling mesh by:\n")
        fout.write("phonopy --qe -t mesh.conf -c pos.in\n")
        fout.write("# Thermal properties can be plotted by:\n")
        fout.write("phonopy --qe -t -p mesh.conf -c pos.in\n")
        fout.write("# calculate Projected DOS and plot it\n")
        fout.write("phonopy --qe -p pdos.conf -c pos.in\n")
        fout.write("# plot band structure\n")
        fout.write("phonopy --qe -p band.conf -c pos.in\n")

    os.system("bash phonopy-analysis.sh")


    import yaml
    with open("band.yaml", 'r') as fin:
        band_yaml = yaml.safe_load(fin)
    
    npath = band_yaml["npath"]
    segment_nqpoint = band_yaml["segment_nqpoint"]
    labels = band_yaml["labels"]
    
    with open("band.data", 'w') as fout:
        #phonon band, different column are different band, but band in different kpath segmentation is not in the same order
        #so we divide it with empty lines, so that there are no weird vertical lines.
        fout.write("#kpoint(segmented) band(1-n)\n")
        nqpoint = band_yaml["nqpoint"]
        nband = band_yaml["phonon"][0]["band"].__len__()
        # actually band_yaml["phonon"].__len__() == nqpoint
        #for i in range(nqpoint):
        #    fout.write("%f" % band_yaml["phonon"][i]["distance"])
        #    for band in band_yaml["phonon"][i]["band"]:
        #        fout.write(" %f" % band["frequency"])
        #    fout.write("\n")

        for s in range(len(segment_nqpoint)):
            if s == 0:
                start = 0
                end = segment_nqpoint[0] - 1                
            else:
                start = sum(segment_nqpoint[:s])
                end = start + segment_nqpoint[s] - 1
            for i in range(start, end + 1):
                fout.write("%f" % band_yaml["phonon"][i]["distance"])
                for band in band_yaml["phonon"][i]["band"]:
                    fout.write(" %f" % band["frequency"])
                fout.write("\n")            
            fout.write("\n")

    locs = []
    labels_for_matplotlib = []
    labels_for_gnuplot = []
    
    labels_for_gnuplot.append(labels[0][0].split("$")[1] if labels[0][0] != "$\Gamma$" else "{/symbol G}")
    labels_for_gnuplot.append(labels[0][1].split("$")[1] if labels[0][1] != "$\Gamma$" else "{/symbol G}")
    labels_for_matplotlib.append(labels[0][0])
    labels_for_matplotlib.append(labels[0][1])
    iqpoint = 0
    locs.append(band_yaml["phonon"][iqpoint]["distance"])
    iqpoint += segment_nqpoint[0] - 1
    locs.append(band_yaml["phonon"][iqpoint]["distance"])
    for ipath in range(1, npath):
        # labels
        if labels[ipath][0] == labels[ipath-1][1]:
            if labels[ipath][1] == "$\Gamma$":
                labels_for_gnuplot.append("{/symbol G}")
            else:
                labels_for_gnuplot.append(labels[ipath][1].split("$")[1])
            labels_for_matplotlib.append(labels[ipath][1])
        else:
            if labels[ipath-1][1] == "$\Gamma$":
                labels_for_gnuplot[-1] = "{/symbole G}" + "|" + labesl[ipath][0].split("$")[1]
            elif labels[ipath][0] == "$\Gamma$":
                labels_for_gnuplot[-1] = labesl[ipath-1][1].split("$")[1] + "|" + "{/symbol G}"
            else:
                labels_for_gnuplot[-1] = labels[ipath-1][1].split("$")[1] + "|" + labels[ipath][0].split("$")[1]
            if labels[ipath][1] == "$\Gamma$":
                labels_for_gnuplot.append("{/symbol G}")
            else:
                labels_for_gnuplot.append(labels[ipath][1].split("$")[1])
            labels_for_matplotlib[-1] = "$" + labels[ipath-1][1].split("$")[1] + "|" + labels[ipath][0].split("$")[1] + "$"
            labels_for_matplotlib.append(labels[ipath][1])
        # locs
        iqpoint += segment_nqpoint[ipath]
        locs.append(band_yaml["phonon"][iqpoint]["distance"])
        
        
    if engine == "matplotlib":
        import numpy as np
        import matplotlib.pyplot as plt
        with open("band.data", 'r') as fin:
            band_data = np.loadtxt(fin)
        # in band.yaml band in different kpath segmentation is not in the same order, we plot each segmentation separately
        # so that there are no weird vertical lines
        for s in range(len(segment_nqpoint)):
            if s == 0:
                start = 0
                end = segment_nqpoint[0] - 1                
            else:
                start = sum(segment_nqpoint[:s])
                end = start + segment_nqpoint[s] - 1            
            for iband in range(nband):
                plt.plot(band_data[start:end+1, 0], band_data[start:end+1, iband+1], color='red', linewidth=1)
        plt.xticks(locs, labels_for_matplotlib)
        plt.xlabel("K")
        plt.ylabel("Frequency (THz)")
        #plt.grid(b=True, which='major')
        #if xrange != None:
        #    plt.xlim(xmin=xrange[0], xmax=xrange[1])
        #if yrange != None:
        #    plt.ylim(ymin=yrange[0], ymax=yrange[1])
        plt.savefig("phonon_band.png")
        plt.close()
    elif engine == "gnuplot":
        with open("band.gnuplot", 'w') as fout:
            fout.write("set terminal gif\n")
            fout.write("set output 'phonon_band.gif'\n")
            fout.write("unset key\n")
            fout.write("set parametric\n")
            #fout.write("set title 'Band Structure (Spin %d)'\n" % (i+1))
            fout.write("set xlabel 'K'\n")
            fout.write("set ylabel 'Frequency(THz)'\n")
            fout.write("set xtics(")
            for j in range(len(labels_for_gnuplot)-1):
                fout.write("'%s' %f, " % (labels_for_gnuplot[j], locs[j]))
            fout.write("'%s' %f)\n" % (labels_for_gnuplot[-1], locs[-1]))
            fout.write("set grid xtics ytics\n")
            fout.write("set autoscale\n")
            #if xrange != None:
            #    fout.write("set xrange [%f:%f]\n" % (xrange[0], xrange[1]))
            #if yrange != None:
            #    fout.write("set yrange [%f:%f]\n" % (yrange[0], yrange[1]))
            fout.write("plot for [i=2:%d:1] 'band.data' using 1:i w l\n" % (nband + 1))
        os.system("gnuplot band.gnuplot")            

    os.chdir("../")
