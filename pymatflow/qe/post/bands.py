#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os

class bands_post:
    """
    """
    def __init__(self, pwxbandsin, bandsxout):
        """
        pwxbandin: 
            the input file for the pw.x band calculation. used to get the special kpoints names.
        bandsxout:
            used to get the x coordinate of the special kpoints

        specialk:
            [{label: "Gamma", coord: [float, float, float], xcoord: float}, ....]
            coord is the kpoint coordinate specified in pw.x band calculation
            xcoord is the correspondin x coordinates for the special kpoint in
            band plot(read from the band.x calculation output: x coordinate)
        """
        with open(pwxbandsin, 'r') as fout:
            self.pwxbandsin = fout.readlines()
        with open(bandsxout, 'r') as fout:
            self.bandsxout = fout.readlines()

        self.specialk = []
        self.bandfile_gnu = None
        self.get_info()

    def get_info(self):
        # get the special kpoint coord and label from pw.x band calculation input file
        nspecialk = 0
        special_k_begin = 0
        special_k_end = 0
        for i in range(len(self.pwxbandsin)):
            if len(self.pwxbandsin[i].split()) == 0:
                continue
            if self.pwxbandsin[i].split()[0] == "K_POINTS":
                nspecialk = int(self.pwxbandsin[i+1].split()[0])
                special_k_begin = i + 2
                special_k_end = i + 1 + nspecialk
        
        for i in range(special_k_begin, special_k_end + 1):
            print(self.pwxbandsin[i].split("#"))
            kpoint = {"label": self.pwxbandsin[i].split("#")[1].split()[0], "coord": [float(self.pwxbandsin[i].split()[0]), float(self.pwxbandsin[i].split()[1]), float(self.pwxbandsin[i].split()[2])], "xcoord": None}
            self.specialk.append(kpoint)

        # get the x coordinate from band.x output file
        xcoord_begin = 0
        for i in range(len(self.bandsxout)):
            if len(self.bandsxout[i].split()) == 0:
                continue
            if self.bandsxout[i].split()[0] == "high-symmetry":
                xcoord_begin = i
                break # break for loop in the first high-symmetry
        for i in range(nspecialk):
            print(self.bandsxout[xcoord_begin+i])
            self.specialk[i]["xcoord"] = float(self.bandsxout[xcoord_begin+i].split()[7])
        #
        # get the xxx.data.gnu file names
        for line in self.bandsxout:
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "Plottable":
                self.bandfile_gnu = line.split()[6] # energy in eV
        #
        # get fermi energy from nscf output
        scfout = "static-scf.out"
        nscfout = "static-nscf.out"
        if os.path.exists(os.path.join("./", nscfout)):
            with open(os.path.join("./", nscfout), 'r') as fin:
                for line in fin:
                    if len(line.split()) == 0:
                        continue
                    if line.split()[0] == "the" and line.split()[1] == "Fermi":
                        self.efermi = float(line.split()[4])
        elif os.path.exists(os.path.join("./", scfout)):
            with open(os.path.join("./", scfout), 'r') as fin:
                for line in fin:
                    if len(line.split()) == 0:
                        continue
                    if line.split()[0] == "the" and line.split()[1] == "Fermi":
                        self.efermi = float(line.split()[4])
        else:
            print("===========================================================\n")
            print("                Warning !!!\n")
            print("===========================================================\n")
            print("BAND structure postprocessing:\n")
            print("must provide nscfout or at least scfout to get Fermi energy\n")
            sys.exit(1)
        # we do not directly shift Efermi to zero in the dataset
        # but only use it to set the gnuplot scripts so that gnuplot
        # script will be responsible for shfiting Efermi to zero.

    def plot_band(self, option="gnuplot"):
        """
        option:
            gnuplot or matplotlib
        """
        if option == "gnuplot":
            with open("bandplot.gp", 'w') as fout:
                fout.write("set terminal gif\n")
                fout.write("set output 'bandstructure.gif'\n")
                fout.write("unset key\n")
                fout.write("set parametric\n")
                
                fout.write("set title 'Bandstructure'\n")
                fout.write("set xlabel 'K'\n")
                fout.write("set ylabel 'Energy(eV)'\n")
                fout.write("set xtics(")
                for point in self.specialk:
                    if point["label"] == "GAMMA":
                        fout.write("'%s' %f, " % ("{/symbol G}", point["xcoord"]))
                    else:
                        fout.write("'%s' %f, " % (point["label"], point["xcoord"]))
                fout.write(")\n")
                fout.write("set grid xtics ytics\n")
                fout.write("set autoscale\n")
                fout.write("# fermi energy shifted to zero by use using 1:($2-%f) in plot function\n" % self.efermi)
                fout.write("# and data in %s file is not modified at all, and is as it is\n" % self.bandfile_gnu)
                fout.write("plot ")
                fout.write("'%s' using 1:($2-%f) w l" % (self.bandfile_gnu, self.efermi))
                # 
                #for i in range(len(self.specialk) - 1):
                #    fout.write(", %f, t" % (self.specialk[i]["xcoord"]))
                #fout.write(", %f, t\n" % self.specialk[-1]["xcoord"])
                #
                fout.write("\n")
            os.system("gnuplot bandplot.gp")

        elif option == "matplotlib":
            pass
