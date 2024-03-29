#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os

"""
About:
    the xxx.bands file generated by siesta when running band calculation
    can be processed by gnubands program(distributred with siesta soruce
    ).
    gnubands extract data from xxx.bands and output an gnuplot format file.
    actually we can easily implement that through python. 

    however, gnubands didn't shift Efermi to zero automatically.
    so the Efermi shifted data file is generated based on the output
    of gnubands.
"""


class BandsPost:
    def __init__(self):
        self.bands_file_lines = None
        self.specialk = None
    
    def process(self, bandsfile="siesta.bands"):
        """
        bandsfile:
            the xxx.bands file generated by siesta
        """
        with open(bandsfile, 'r') as fin:
            self.bands_file_lines = fin.readlines()

        self.get_specialk()
            
        # use gnubands to convert xxx.bands to gnuplot applicable data
        os.system("gnubands %s > %s.gnuplot.data" % (bandsfile, bandsfile))

        # shift efermi to 0
        efermi = float(self.bands_file_lines[0].split()[0])
        with open("%s.gnuplot.data" % bandsfile, 'r') as fin:
            lines = fin.readlines()
        with open("%s.gnuplot.shifted.data" % bandsfile, 'w') as fout:
            fout.write("# ===================================================\n")
            fout.write("# Pymatflow:\n")
            fout.write("# based on gnubands, but E fermi is shfited to zero!\n")
            fout.write("# ===================================================\n")
            for line in lines:
                if len(line.split()) == 0:
                    fout.write(line)
                    continue
                if line.split()[0] == "#":
                    fout.write(line)
                else:
                    fout.write("\t%.6f\t%.6f\t%d\n" % (float(line.split()[0]), float(line.split()[1]) - efermi, 1))
        #

    def get_specialk(self):
        """
        self.specialk:
            [{"label": "GAMMA", "xcoord": float}, ...]
        Note:
            extract the x coordinate and label of special kpoint
        """
        nspecialk = 0
        special_k_begin = 0
        for i in range(len(self.bands_file_lines)):
            if len(self.bands_file_lines[i].split()) == 1 and i != 0:
                nspecialk = int(self.bands_file_lines[i].split()[0])
                special_k_begin = i + 1
        self.specialk = []
        for i in range(special_k_begin, special_k_begin+nspecialk):
            self.specialk.append({"label": self.bands_file_lines[i].split()[1], "xcoord": float(self.bands_file_lines[i].split()[0])})

    def plot_bands(self, bandsfile="siesta.bands", option="gnuplot"):
        """
        option:
            use gnuplot or matplotlib to plot the band structure
        """
        if option == "gnuplot":
            with open("bandplot.gp", 'w') as fout:
                fout.write("set terminal gif\n")
                fout.write("set output 'bandstructure.gif'\n")
                fout.write("set title 'Bandstructure'\n")
                fout.write("set xlabel 'K'\n")
                fout.write("set ylabel 'Energy(eV)'\n")
                fout.write("unset key\n")
                fout.write("set xtics(")
                for point in self.specialk:
                    if point["label"] == 'GAMMA':
                        fout.write("%s %f, " % ("{/symbol G}", point["xcoord"]))
                    else:
                        fout.write("%s %f, " % (point["label"], point["xcoord"])) # minus 1, because in gnuplot x start with 0
                fout.write(")\n")
                fout.write("plot '%s.gnuplot.shifted.data' u 1:2  w l \n" % (bandsfile))
            os.system("gnuplot bandplot.gp")

        elif option == "matplotlib":
            pass
    def export(self, directory, option="gnuplot"):
        os.chdir(directory)
        os.system("mkdir -p post-processing")
        os.chdir("post-processing")
        self.plot_bands(bandsfile="../siesta.bands", option=option)
        os.chdir("../../")
