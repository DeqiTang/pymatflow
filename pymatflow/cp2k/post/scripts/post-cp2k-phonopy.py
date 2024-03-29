#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.base.xyz import base_xyz

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-cp2k-phonopy",
            help="directory of phonopy running")

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
    
    xyz = base_xyz()
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

    #
    # get the disps information
    os.chdir(args.directory)
    os.system("mkdir -p post-processing")

    os.system("ls | grep 'phonon-supercell-' > ./post-processing/geo.data")
    disps = []
    with open("post-processing/geo.data", 'r') as fin:
        for line in fin:
            disps.append(line.split(".")[0].split("-")[2])


    # generate the result analysis bash script and necessary config files

    with open("post-processing/mesh.conf", 'w') as fout:
        fout.write("ATOM_NAME =")
        for element in xyz.specie_labels:
            fout.write(" %s" % element)
        fout.write("\n")
        fout.write("DIM = %d %d %d\n" % (args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))
        fout.write("MP = %d %d %d\n" % (args.mp[0], args.mp[1], args.mp[2]))



    with open("post-processing/pdos.conf", 'w') as fout:
        fout.write("ATOM_NAME =")
        for element in xyz.specie_labels:
            fout.write(" %s" % element)
        fout.write("\n")
        fout.write("DIM = %d %d %d\n" % (args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))
        fout.write("MP = %d %d %d\n" % (args.mp[0], args.mp[1], args.mp[2]))
        fout.write("PDOS = 1 2, 3 4 5 5\n")

    with open("post-processing/band.conf", 'w') as fout:
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
        for qpoint in qpath:
            if qpoint[4] == None:
                fout.write(" %f %f %f" % (qpoint[0], qpoint[1], qpoint[2]))
            elif qpoint[4] == "|":
                fout.write(" %f %f %f," % (qpoint[0], qpoint[1], qpoint[2]))
            else:
                pass
        fout.write("\n")
        fout.write("BAND_LABELS =")
        for qpoint in qpath:
            if qpoint[4] == None:
                if qpoint[3].upper() == "GAMMA":
                    fout.write(" $\Gamma$")
                else:
                    fout.write(" $%s$" % qpoint[3])
            elif qpoint[4] == "|":
                if qpoint[3].upper() == "GAMMA":
                    fout.write(" $\Gamma$,")
                else:
                    fout.write(" $%s$," % qpoint[3])
            else:
                pass
        fout.write("\n")

    with open("post-processing/phonopy-analysis.sh", 'w') as fout:
        inp_name = "phonon.inp"
        fout.write("#!/bin/bash\n\n")
        fout.write("cp ../phonopy_disp.yaml ./\n")
        fout.write("# get the FORCE_SETS\n")
        base_project_name = "ab-initio"
        fout.write("phonopy --cp2k -f ../%s-supercell-{001..%s}-forces-1_0.xyz\n" % (base_project_name, disps[-1]))
        fout.write("# plot The density of states (DOS)\n")
        fout.write("phonopy --cp2k -p mesh.conf -c ../%s\n" % inp_name)
        fout.write("# Thermal properties are calculated with the sampling mesh by:\n")
        fout.write("phonopy --cp2k -t mesh.conf -c ../%s\n" % inp_name)
        fout.write("# Thermal properties can be plotted by:\n")
        fout.write("phonopy --cp2k -t -p mesh.conf -c ../%s\n" % inp_name)
        fout.write("# calculate Projected DOS and plot it\n")
        fout.write("phonopy --cp2k -p pdos.conf -c ../%s\n" % inp_name)
        fout.write("# get the band structure\n")
        fout.write("phonopy --cp2k -c ../%s -p band.conf\n" % inp_name)

    os.system("cd post-processing; bash phonopy-analysis.sh; cd ../")

    import yaml
    with open("post-processing/band.yaml", 'r') as fin:
        band_yaml = yaml.safe_load(fin)
    
    npath = band_yaml["npath"]
    segment_nqpoint = band_yaml["segment_nqpoint"]
    labels = band_yaml["labels"]
    
    with open("post-processingband.data", 'w') as fout:
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
        with open("post-processing/band.data", 'r') as fin:
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
        plt.savefig("post-processing/phonon_band.png")
        plt.close()
    elif engine == "gnuplot":
        with open("post-processing/band.gnuplot", 'w') as fout:
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
        os.system("cd post-processing; gnuplot band.gnuplot; cd ../")            
        
    os.chdir("../")
