"""
post_bands:
    post_bands extract data from static-o_DS3_EBANDS.agr and it will build
    the kpoints length: xcoord_k from the high symmetry line and the corresponding
    basis for reciprocal space.

    b1 = 1 / a1, b2 = 1 / a2 and b3 = 1 / a3.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

class post_bands:
    def __init__(self):
        pass


    def get_xcoord_k(self, kpath, cell):
        """
        Note:
            xcoord_k is the x axis of the band structure plot
            let's see how it is build from kpoints and the
            crystal lattice or reciprocal lattice.
        """
        self.kpath = kpath
        self.xcoord_k = []
        b1 = 1 / np.sqrt(cell[0][0]**2 + cell[0][1]**2 + cell[0][2]**2)
        b2 = 1 / np.sqrt(cell[1][0]**2 + cell[1][1]**2 + cell[1][2]**2)
        b3 = 1 / np.sqrt(cell[2][0]**2 + cell[2][1]**2 + cell[2][2]**2)
        # actually you will find that in vasp b1=1/a1, b2=1/a2, b3=1/a3

        self.xcoord_k.append(0.0000000)
        for i in range(len(self.kpath) - 1):
            # the step in the xcoord_k for each segment is different and it is actually
            # the distance between the two high symmetry kpoint in unit of reciprocal coordinates
            # divided by the conneciting number kpath[i][4]
            if self.kpath[i][4] != "|":
                step = 0
                delta_b_1 = b1*(self.kpath[i+1][0] - self.kpath[i][0])
                delta_b_2 = b2*(self.kpath[i+1][1] - self.kpath[i][1])
                delta_b_3 = b3*(self.kpath[i+1][2] - self.kpath[i][2])
                step = np.sqrt(delta_b_1**2+delta_b_2**2+delta_b_3**2) / (self.kpath[i][4])

                for j in range(self.kpath[i][4]):
                    self.xcoord_k.append(self.xcoord_k[-1]+step)
            else:
                self.xcoord_k.append(self.xcoord_k[-1])
        # label for plot
        self.locs = []
        self.labels_for_matplotlib = []
        self.labels_for_gnuplot = []

        self.locs.append(0.0000000)
        nk = 0
        print("%d\n" % nk)
        for i in range(len(self.kpath)-1):
            if self.kpath[i][4] != "|":
                nk = nk + self.kpath[i][4]
                self.locs.append(self.xcoord_k[nk])
                print("%d\n" % nk)
            else:
                nk = nk + 1

        self.labels_for_matplotlib.append(r"$%s$" % self.kpath[0][3].upper() if self.kpath[0][3].upper() != "GAMMA" else r"$\Gamma$")
        self.labels_for_gnuplot.append("%s" % self.kpath[0][3].upper() if self.kpath[0][3].upper() != "GAMMA" else "{/symbol G}")
        for i in range(1, len(self.kpath)):
            if self.kpath[i-1][4] != "|":
                self.labels_for_matplotlib.append(r"$%s$" % self.kpath[i][3].upper() if self.kpath[i][3].upper() != "GAMMA" else r"$\Gamma$")
                self.labels_for_gnuplot.append("%s" % self.kpath[i][3].upper() if self.kpath[i][3].upper() != "GAMMA" else "{/symbol G}")
            else:
                self.labels_for_matplotlib[-1] = r"$%s | %s$" % (self.labels_for_matplotlib[-1].split("$")[1], self.kpath[i][3].upper())
                self.labels_for_gnuplot[-1] = "%s | %s" % (self.labels_for_gnuplot[-1], self.kpath[i][3].upper())


    def get_ebands_agr(self, filepath="static-o_DS3_EBANDS.agr"):
        with open(filepath, 'r') as fin:
            self.lines = fin.readlines()
        # get the band energy
        # in xxx_EBANDS.agr, energy are in unit of eV, and fermi energy are already shfited to 0

        # first check the magnetic_status
        for line in self.lines:
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "#" and line.split()[1] == "mband:":
                self.mband = int(line.split()[2].split(",")[0])
                self.nkpt = int(line.split()[4].split(",")[0])
                self.nsppol = int(line.split()[6].split(",")[0])
                self.nspinor = int(line.split()[8])
        # get the eigenval (in agr, efermi is shfited to 0 already)
        self.energies_agr = []
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0] == "@type" and self.lines[i].split()[1].split("\n")[0] == "xy":
                band = []
                for j in range(self.nkpt):
                    band.append(float(self.lines[i+j+1].split()[1]))
                self.energies_agr.append(band)

    def _plot_band_matplotlib(self, bandrange=[0, 1.0]):
        """
        :param bandrange:
            a list of two values(between 0 and 1) defining the percentage
            of bands to plot.
            plotrange[0]: left boundary of the nth band to plot
            plotrange[1]: right boundary of the nth band to plot
            default is plotrange[0] = 0, plotrange[1], in which case
            all the band available will be plot.
            Be aware that the range if not for energy but for band number
        :param imagebase: image file name(not full)
        """

        if self.nsppol == 1:
            band_min = int(bandrange[0] * self.mband)
            band_max = int(bandrange[1] * self.mband)
            for i in range(band_min, band_max, 1):
                plt.plot(self.xcoord_k, self.energies_agr[i], color='blue', linewidth=1)
            plt.xticks(self.locs, self.labels_for_matplotlib)
            plt.title("Band Structure")
            plt.xlabel("K")
            plt.ylabel("Energy(eV)")
            plt.grid(b=True, which='major')
            plt.savefig("band-structure-spin-unpolarized.png")

        if self.nsppol == 2:
            # half of self.energies_agr are spin up, and half are spin down
            band_min = int(bandrange[0] * self.mband)
            band_max = int(bandrange[1] * self.mband)
            # spin up
            for i in range(band_min, band_max, 1):
                plt.plot(self.xcoord_k, self.energies_agr[i])
            plt.title("Band Structure(Spin Up)")
            plt.xlabel("K")
            plt.ylabel("Energy(eV)")
            plt.xticks(self.locs, self.labels_for_matplotlib)
            plt.grid(b=True, which='major')
            plt.savefig("band-structure-spin-polarized-1.png")
            plt.close()
            # spin down
            for i in range(int(band_min+self.mband), int(band_max+self.mband), 1):
                plt.plot(self.xcoord_k, self.energies_agr[i])
            plt.title("Band Structure(Spin Down)")
            plt.xlabel("K")
            plt.ylabel("Energy(eV)")
            plt.xticks(self.locs, self.labels_for_matplotlib)
            plt.grid(b=True, which='major')
            plt.savefig("band-structure-spin-polarized-2.png")
            plt.close()
            # all in one picture
            for i in range(band_min, band_max, 1):
                plt.plot(self.xcoord_k, self.energies_agr[i], color="blue", linewidth=1)
            for i in range(int(band_min+self.mband), int(band_max+self.mband), 1):
                plt.plot(self.xcoord_k, self.energies_agr[i], color="red", linewidth=1)
            plt.title("Band Structure(Spin Up&Down)")
            plt.xlabel("K")
            plt.ylabel("Energy(eV)")
            plt.xticks(self.locs, self.labels_for_matplotlib)
            plt.grid(b=True, which='major')
            plt.savefig("band-structure-spin-polarized-all.png")
            plt.close()


    def _plot_band_gnuplot(self, bandrange=[0, 1.0]):
        """
        :param bandrange:
            a list of two values(between 0 and 1) defining the percentage
            of bands to plot.
            plotrange[0]: left boundary of the nth band to plot
            plotrange[1]: right boundary of the nth band to plot
            default is plotrange[0] = 0, plotrange[1], in which case
            all the band available will be plot.
            Be aware that the range if not for energy but for band number
        :param imagebase: image file name(not full)
        """

        if self.nsppol == 1:
            band_min = int(bandrange[0] * self.mband)
            band_max = int(bandrange[1] * self.mband)
            with open("all-bands-spin-unpolarized.data", 'w') as fout:
                fout.write("# band structure extracted from xxx_EBANDS.agr\n")
                fout.write("# efermi shfited to 0 already\n")
                for i in range(self.mband):
                    for j in range(len(self.xcoord_k)):
                        fout.write("%f %f\n" % (self.xcoord_k[j], self.energies_agr[i][j]))
                    fout.write("\n")
            with open("specified-bands-spin-unpolarized.data", 'w') as fout:
                fout.write("# band structure extracted from ***_EBANDS.agr\n")
                fout.write("# efermi shifted to 0 already\n")
                for i in range(band_min, band_max, 1):
                    for j in range(len(self.xcoord_k)):
                        fout.write("%f %f\n" % (self.xcoord_k[j], self.energies_agr[i][j]))
                    fout.write("\n")

            with open("all-bands-spin-unpolarized.gnuplot", 'w') as fout:
                fout.write("set terminal gif\n")
                fout.write("set output 'all-bands-spin-unpolarized.gif'\n")
                fout.write("unset key\n")
                fout.write("set parametric\n")
                fout.write("set title 'Band Structure'\n")
                fout.write("set xlabel 'K'\n")
                fout.write("set ylabel 'Energy(eV)'\n")
                fout.write("set xtics(")
                for i in range(len(self.labels_for_gnuplot)-1):
                    fout.write("'%s' %f, " % (self.labels_for_gnuplot[i], self.locs[i]))
                fout.write("'%s' %f)\n" % (self.labels_for_gnuplot[-1], self.locs[-1]))
                fout.write("set grid xtics ytics\n")
                fout.write("set autoscale\n")
                fout.write("plot 'all-bands-spin-unpolarized.data' using 1:2 w l\n")
            os.system("gnuplot all-bands-spin-unpolarized.gnuplot")

            with open("specified-bands-spin-unpolarized.gnuplot", 'w') as fout:
                fout.write("set terminal gif\n")
                fout.write("set output 'specified-bands-spin-unpolarized.gif'\n")
                fout.write("unset key\n")
                fout.write("set parametric\n")
                fout.write("set title 'Band Structure'\n")
                fout.write("set xlabel 'K'\n")
                fout.write("set ylabel 'Energy(eV)'\n")
                fout.write("set xtics(")
                for i in range(len(self.labels_for_gnuplot)-1):
                    fout.write("'%s' %f, " % (self.labels_for_gnuplot[i], self.locs[i]))
                fout.write("'%s' %f)\n" % (self.labels_for_gnuplot[-1], self.locs[-1]))
                fout.write("set grid xtics ytics\n")
                fout.write("set autoscale\n")
                fout.write("plot 'specified-bands-spin-unpolarized.data' using 1:2 w l\n")
            os.system("gnuplot specified-bands-spin-unpolarized.gnuplot")

        if self.nsppol == 2:
            band_min = int(bandrange[0] * self.mband)
            band_max = int(bandrange[1] * self.mband)

            with open("all-bands-spin-polarized-spin-1.data", 'w') as fout:
                fout.write("# band structure extracted from xxx_EBANDS.agr\n")
                fout.write("# efermi shfited to 0 already\n")
                for i in range(self.mband):
                    for j in range(len(self.xcoord_k)):
                        fout.write("%f %f\n" % (self.xcoord_k[j], self.energies_agr[i][j]))
                    fout.write("\n")

            with open("all-bands-spin-polarized-spin-2.data", 'w') as fout:
                fout.write("# band structure extracted from xxx_EBANDS.agr\n")
                fout.write("# efermi shifted to 0 already\n")
                for i in range(self.mband):
                    for j in range(len(self.xcoord_k)):
                        fout.write("%f %f\n" % (self.xcoord_k[j], self.energies_agr[self.mband+i][j]))
                    fout.write("\n")

            with open("specified-bands-spin-polarized-spin-1.data", 'w') as fout:
                fout.write("# band structure extracted from xxx_EBANDS.agr\n")
                fout.write("# efermi shifted to 0 already\n")
                for i in range(band_min, band_max, 1):
                    for j in range(len(self.xcoord_k)):
                        fout.write("%f %f\n" % (self.xcoord_k[j], self.energies_agr[i][j]))
                    fout.write("\n")

            with open("specified-bands-spin-polarized-spin-2.data", 'w') as fout:
                fout.write("# band structure extracted from xxx_EBANDS.agr\n")
                fout.write("# efermi shifted to 0 already\n")
                for i in range(band_min, band_max, 1):
                    for j in range(len(self.xcoord_k)):
                        fout.write("%f %f\n" % (self.xcoord_k[j], self.energies_agr[self.mband+i][j]))
                    fout.write("\n")



            with open("all-bands-spin-polarized-spin-1.gnuplot", 'w') as fout:
                fout.write("set terminal gif\n")
                fout.write("set output 'all-bands-spin-polarized-spin-1.gif'\n")
                fout.write("unset key\n")
                fout.write("set parametric\n")
                fout.write("set title 'Band Structure'\n")
                fout.write("set xlabel 'K'\n")
                fout.write("set ylabel 'Energy(eV)'\n")
                fout.write("set xtics(")
                for i in range(len(self.labels_for_gnuplot)-1):
                    fout.write("'%s' %f, " % (self.labels_for_gnuplot[i], self.locs[i]))
                fout.write("'%s' %f)\n" % (self.labels_for_gnuplot[-1], self.locs[-1]))
                fout.write("set grid xtics ytics\n")
                fout.write("set autoscale\n")
                fout.write("plot 'all-bands-spin-polarized-spin-1.data' using 1:2 w l\n")
            os.system("gnuplot all-bands-spin-polarized-spin-1.gnuplot")


            with open("all-bands-spin-polarized-spin-2.gnuplot", 'w') as fout:
                fout.write("set terminal gif\n")
                fout.write("set output 'all-bands-spin-polarized-spin-2.gif'\n")
                fout.write("unset key\n")
                fout.write("set parametric\n")
                fout.write("set title 'Band Structure'\n")
                fout.write("set xlabel 'K'\n")
                fout.write("set ylabel 'Energy(eV)'\n")
                fout.write("set xtics(")
                for i in range(len(self.labels_for_gnuplot)-1):
                    fout.write("'%s' %f, " % (self.labels_for_gnuplot[i], self.locs[i]))
                fout.write("'%s' %f)\n" % (self.labels_for_gnuplot[-1], self.locs[-1]))
                fout.write("set grid xtics ytics\n")
                fout.write("set autoscale\n")
                fout.write("plot 'all-bands-spin-polarized-spin-2.data' using 1:2 w l\n")
            os.system("gnuplot all-bands-spin-polarized-spin-2.gnuplot")


            with open("specified-bands-spin-polarized-spin-1.gnuplot", 'w') as fout:
                fout.write("set terminal gif\n")
                fout.write("set output 'specified-bands-spin-polarized-spin-1.gif'\n")
                fout.write("unset key\n")
                fout.write("set parametric\n")
                fout.write("set title 'Band Structure'\n")
                fout.write("set xlabel 'K'\n")
                fout.write("set ylabel 'Energy(eV)'\n")
                fout.write("set xtics(")
                for i in range(len(self.labels_for_gnuplot)-1):
                    fout.write("'%s' %f, " % (self.labels_for_gnuplot[i], self.locs[i]))
                fout.write("'%s' %f)\n" % (self.labels_for_gnuplot[-1], self.locs[-1]))
                fout.write("set grid xtics ytics\n")
                fout.write("set autoscale\n")
                fout.write("plot 'specified-bands-spin-polarized-spin-1.data' using 1:2 w l\n")
            os.system("gnuplot specified-bands-spin-polarized-spin-1.gnuplot")


            with open("specified-bands-spin-polarized-spin-2.gnuplot", 'w') as fout:
                fout.write("set terminal gif\n")
                fout.write("set output 'specified-bands-spin-polarized-spin-2.gif'\n")
                fout.write("unset key\n")
                fout.write("set parametric\n")
                fout.write("set title 'Band Structure'\n")
                fout.write("set xlabel 'K'\n")
                fout.write("set ylabel 'Energy(eV)'\n")
                fout.write("set xtics(")
                for i in range(len(self.labels_for_gnuplot)-1):
                    fout.write("'%s' %f, " % (self.labels_for_gnuplot[i], self.locs[i]))
                fout.write("'%s' %f)\n" % (self.labels_for_gnuplot[-1], self.locs[-1]))
                fout.write("set grid xtics ytics\n")
                fout.write("set autoscale\n")
                fout.write("plot 'specified-bands-spin-polarized-spin-2.data' using 1:2 w l\n")
            os.system("gnuplot specified-bands-spin-polarized-spin-2.gnuplot")


    def plot_band(self, option="matplotlib", bandrange=[0, 1.0]):
        """
        :parama option:
            gnuplot or matplotlib
        :param bandrange:
            a list of two values(between 0 and 1) defining the percentage
            of bands to plot.
            plotrange[0]: left boundary of the nth band to plot
            plotrange[1]: right boundary of the nth band to plot
            default is plotrange[0] = 0, plotrange[1], in which case
            all the band available will be plot.
            Be aware that the range if not for energy but for band number
        :param imagebase: image file name(not full)
        """
        if option == "matplotlib":
            self._plot_band_matplotlib(bandrange=bandrange)
        elif option == "gnuplot":
            self._plot_band_gnuplot(bandrange=bandrange)
        #

    def export(self, directory="tmp-abinit-static", bandrange=[0, 1], option="matplotlib"):
        """
        :parama option:
            gnuplot or matplotlib
        :param bandrange:
            a list of two values(between 0 and 1) defining the percentage
            of bands to plot.
            plotrange[0]: left boundary of the nth band to plot
            plotrange[1]: right boundary of the nth band to plot
            default is plotrange[0] = 0, plotrange[1], in which case
            all the band available will be plot.
            Be aware that the range if not for energy but for band number
        """
        os.system("mkdir -p %s/post-processing" % directory)
        os.chdir(os.path.join(directory, "post-processing"))
        self.plot_band(option=option,  bandrange=bandrange)
        os.chdir("../../")
