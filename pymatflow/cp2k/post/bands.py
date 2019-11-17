#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os

class bands_post:
    """
    """
    def __init__(self):
        """
        """
        self.get_info()

    def get_info(self):
        pass

    def plot_band(self, bandsfile="bands.bs", option="gnuplot"):
        """
        option:
            gnuplot or matplotlib
        specialk:
            [{"label": "GAMMA", "coord": [float, float, float], "k-number"}, ...]
            where k-number means the speical k is the n-th kpoint in the total kpoints
        """
        specialk = []
        with open(bandsfile, 'r') as fin:
            lines = fin.readlines()
        for line in lines:
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "POINT":
                specialk.append({"label": line.split()[2], "coord": [float(line.split()[3]), float(line.split()[4]), float(line.split()[5])], "k-number": None})
        # get the k-number for each special kpoint
        for line in lines:
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "Nr.":
                for kpoint in specialk:
                    if [float(line.split()[5]), float(line.split()[6]), float(line.split()[7])] == kpoint["coord"] and kpoint["k-number"] is None:
                        kpoint["k-number"] = int(line.split()[1])
                        break
                        # must break here to avoid set the next special kpoint with the same label and coord
        #
        # get the number of bands
        n_bands = 0
        for i in range(len(lines)):
            if len(lines[i].split()) == 0:
                continue
            if lines[i].split()[0] == "Nr.":
                n_bands = int(lines[i+1].split()[0])

        #
        if option == "gnuplot":
            os.system("cp2k_bs2csv.py %s" % bandsfile) 
            # cp2k_bs2csv.py is provided by cp2k.org and I put it in this library too.
            # see https://www.cp2k.org/exercises:2018_uzh_cmest:pdos
            with open("bandplot.gp", 'w') as fout:
                fout.write("set terminal gif\n")
                fout.write("set output 'bandstructure.gif'\n")
                fout.write("set title 'Bandstructure'\n")
                #fout.write("set xlabel 'Kpoint number'\n") # equal to nmber of kpoints
                fout.write("set ylabel 'Energy'\n")
                fout.write("unset key\n")
                fout.write("set xtics(")
                for point in specialk:
                    fout.write("'%s' %d, " % (point["label"], point["k-number"]-1)) # minus 1, because in gnuplot x start with 0
                fout.write(")\n")
                fout.write("plot for [i=4:%d] '%s.set-1.csv' u 0:i w l \n" % (n_bands+3, bandsfile))
            os.system("gnuplot bandplot.gp")

        elif option == "matplotlib":
            pass
