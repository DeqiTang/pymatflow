#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import numpy as np
import matplotlib.pyplot as plt

class BandsPost:
    """
    """
    def __init__(self):
        """
        efermi in static-scf.out is in unit of a.u.
        but eigenvalue in bands.bs is in unit of eV.
        so we convert efermi from a.u. to eV
        """
        pass
        self.ha_to_ev = 27.211324570273 # needed to convert Fermi from [a.u.] to eV

    def get_efermi(self, static_out="static-scf.out"):
        with open(static_out, 'r') as fin:
            lines = fin.readlines()
            
        self.efermi = None
        for line in lines:
            if "Fermi energy:" in line:
                self.efermi = float(line.split()[2]) * self.ha_to_ev # now become eV
        return self.efermi

    def get_kpath_and_bands(self, kpath, cell, bands="bands.bs"):
        # bands.bs
        with open(bands, 'r') as fin:
            self.bands_bs = fin.readlines()

        # kpath
        self.kpath = kpath

        self.get_xcoord_k(cell=cell)
        self.get_eigenval()

        self.locs = []
        self.labels_for_matplotlib = []
        self.labels_for_gnuplot = []

        self.locs.append(self.xcoord_k[0])
        for i in range(1, len(self.xcoord_k) - 1):
            if self.xcoord_k[i] == self.xcoord_k[i-1]:
                self.locs.append(self.xcoord_k[i])
        self.locs.append(self.xcoord_k[-1])

        self.labels_for_matplotlib.append(r"$%s$" % self.kpath[0][3].upper() if self.kpath[0][3].upper() != "GAMMA" else r"$\Gamma$")
        self.labels_for_gnuplot.append("%s" % self.kpath[0][3].upper() if self.kpath[0][3].upper() != "GAMMA" else "{/symbol G}")
        for i in range(1, len(self.kpath)):
            if self.kpath[i-1][4] != "|":
                self.labels_for_matplotlib.append(r"$%s$" % self.kpath[i][3].upper() if self.kpath[i][3].upper() != "GAMMA" else r"$\Gamma$")
                self.labels_for_gnuplot.append("%s" % self.kpath[i][3].upper() if self.kpath[i][3].upper() != "GAMMA" else "{/symbol G}")
            else:
                self.labels_for_matplotlib[-1] = r"$%s | %s$" % (self.labels_for_matplotlib[-1].split("$")[1], self.kpath[i][3].upper())
                self.labels_for_gnuplot[-1] = "%s | %s" % (self.labels_for_gnuplot[-1], self.kpath[i][3].upper())

        #print(self.locs)
        #print(self.labels_for_matplotlib)
        #print(self.labels_for_gnuplot)
        self.get_eigenval()


    def get_xcoord_k(self, cell):
        """
        Note:
            xcoord_k is the x axis of the band structure plot
            let's see how it is build from kpoints and the
            crystal lattice or reciprocal lattice.
        """
        self.xcoord_k = []
        # get the lattice parameter and the reciprocal basis
        a1 = np.sqrt(cell[0][0]**2 + cell[0][1]**2 + cell[0][2]**2)

        a2 = np.sqrt(cell[1][0]**2 + cell[1][1]**2 + cell[1][2]**2)
        
        a3 = np.sqrt(cell[2][0]**2 + cell[2][1]**2 + cell[2][2]**2)

        b1 = 1 / a1

        b2 = 1 / a2

        b3 = 1 / a3

        # actually you will find that in vasp b1=1/a1, b2=1/a2, b3=1/a3
        # now we use the reciprocal lattice constant and the kpoints in crystal coordinate to build the xcoord_k


        # in the past, we read kpoint from vasprun.xml but now we directly read kpoint from self.kpath!
        self.xcoord_k.append(0.0000000)
        for i in range(1, len(self.kpath)):
            if self.kpath[i-1][4] != "|":
                step = 0
                delta_b_1 = b1 * (self.kpath[i][0] - self.kpath[i-1][0])
                delta_b_2 = b2 * (self.kpath[i][1] - self.kpath[i-1][1])
                delta_b_3 = b3 * (self.kpath[i][2] - self.kpath[i-1][2])
                step = np.sqrt(delta_b_1**2+delta_b_2**2+delta_b_3**2) / (self.kpath[i-1][4] - 1)

                for j in range(self.kpath[i-1][4]-1):
                    self.xcoord_k.append(self.xcoord_k[-1]+step)
                if i == len(self.kpath) - 1:
                    pass
                else:
                    self.xcoord_k.append(self.xcoord_k[-1])

    def get_eigenval(self):
        """
        eigenval in xxx.bs file is in unit of eV by default so we need not to convert it
        """
        # first check the magnetic_status
        self.magnetic_status = "spin-unpolarized"
        for line in self.bands_bs:
            if "Spin 2" in line:
                self.magnetic_status = "spin-polarized"
                break

        # get the eigenval
        """
        in case of spin-unpolarized
        self.eigenval = {
            "spin_1": [
                {
                    "energy": [],
                },
                {
                    "energy": [],
                }
                ...... 
            ], # length equals to number of kpoints
        }

        in case of spin-polarized
        self.eigenval = {
            "spin_1": [
                {
                    "energy": [],
                },
                {
                    "energy": [],
                }
                ...... 
            ], # length equals to number of kpoints
            "spin_2": [
                {
                    "energy": [],
                },
                {
                    "energy": [],
                }
                ...... 
            ], # length equals to number of kpoints        
        }        
        """
        
        if self.magnetic_status == "spin-unpolarized":
            self.eigenval = {"spin_1": []}
            for i in range(len(self.bands_bs)):
                if "Nr." in self.bands_bs[i]:
                    eigen = {"energy": []}
                    nband = int(self.bands_bs[i+1].split()[0])
                    nline = np.ceil(nband / 4)
                    ncol_last_line = nband % 4
                    j = 2
                    while j <= (nline+1):
                        if j != (nline+1):
                            for k in range(4):
                                eigen["energy"].append(float(self.bands_bs[i+j].split()[k]))
                        else:
                            for k in range(ncol_last_line):
                                eigen["energy"].append(float(self.bands_bs[i+j].split()[k]))
                        j += 1
                    self.eigenval["spin_1"].append(eigen)
        elif self.magnetic_status == "spin-polarized":
            self.eigenval = {"spin_1": [], "spin_2": []}
            for i in range(len(self.bands_bs)):
                if "Nr." in self.bands_bs[i] and "Spin 1" in self.bands_bs[i]:
                    eigen = {"energy": []}
                    nband = int(self.bands_bs[i+1].split()[0])
                    nline = np.ceil(nband / 4)
                    ncol_last_line = nband % 4
                    j = 2
                    while j <= (nline+1):
                        if j != (nline+1):
                            for k in range(4):
                                eigen["energy"].append(float(self.bands_bs[i+j].split()[k]))
                        else:
                            for k in range(ncol_last_line):
                                eigen["energy"].append(float(self.bands_bs[i+j].split()[k]))
                        j += 1
                    self.eigenval["spin_1"].append(eigen)
                if "Nr." in self.bands_bs[i] and "Spin 2" in self.bands_bs[i]:
                    eigen = {"energy": []}
                    nband = int(self.bands_bs[i+1].split()[0])
                    nline = np.ceil(nband / 4)
                    ncol_last_line = nband % 4
                    j = 2
                    while j <= (nline+1):
                        if j != (nline+1):
                            for k in range(4):
                                eigen["energy"].append(float(self.bands_bs[i+j].split()[k]))
                        else:
                            for k in range(ncol_last_line):
                                eigen["energy"].append(float(self.bands_bs[i+j].split()[k]))
                        j += 1
                    self.eigenval["spin_2"].append(eigen)
                                        
    def _plot_band_matplotlib(self, bandrange=[0, 1.0], xrange=None, yrange=None):
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

        if self.magnetic_status == "spin-unpolarized":
            spin_n = 1
        elif self.magnetic_status == "spin-polarized":
            spin_n = 2

        nband = len(self.eigenval["spin_1"][0]["energy"])
        band_min = int(bandrange[0] * nband)
        band_max = int(bandrange[1] * nband)

        # export the data in gnuplot format
        for i in range(spin_n):
            with open("all-bands-%s-spin-%d.data" % (self.magnetic_status, i+1), 'w') as fout:
                fout.write("# band structure extracted from bands.bs\n")
                if self.efermi == None:
                    fout.write("# efermi: smearing is not used so fermi energy is not calculated in cp2k\n")
                else:
                    fout.write("# efermi: %f (eV) and eigenvalues are in unit of eV\n" % self.efermi)
                for k in range(nband):
                    for j in range(len(self.xcoord_k)):
                        fout.write("%f %f\n" % (self.xcoord_k[j], self.eigenval["spin_%d" % (i+1)][j]["energy"][k]))
                    fout.write("\n")
            with open("specified-bands-%s-spin-%d.data" % (self.magnetic_status, i+1), 'w') as fout:
                fout.write("# band structure extracted from bands.bs\n")
                if self.efermi == None:
                    fout.write("# efermi: smearing is not used so fermi energy is not calculated in cp2k\n")
                else:
                    fout.write("# efermi: %f (eV) and eigenvalues are in unit of eV\n" % self.efermi)
                for k in range(band_min, band_max, 1):
                    for j in range(len(self.xcoord_k)):
                        fout.write("%f %f\n" % (self.xcoord_k[j], self.eigenval["spin_%d" % (i+1)][j]["energy"][k]))
                    fout.write("\n")
        # make the plot
        for i in range(spin_n):
            for j in range(band_min, band_max, 1):
                if self.efermi == None:
                    plt.plot(self.xcoord_k, [self.eigenval["spin_%d" % (i+1)][k]["energy"][j] for k in range(len(self.eigenval["spin_%d" % (i+1)]))], color='blue', linewidth=1)
                else:
                    plt.plot(self.xcoord_k, [self.eigenval["spin_%d" % (i+1)][k]["energy"][j] - self.efermi for k in range(len(self.eigenval["spin_%d" % (i+1)]))], color='blue', linewidth=1)
            plt.xticks(self.locs, self.labels_for_matplotlib)
            #plt.title("Band Structure (Spin %d)" % (i+1))
            plt.xlabel("K")
            plt.ylabel("Energy(eV)")
            plt.grid(b=True, which='major')
            if xrange != None:
                plt.xlim(xmin=xrange[0], xmax=xrange[1])
            if yrange != None:
                plt.ylim(ymin=yrange[0], ymax=yrange[1])
            plt.savefig("band-structure-%s-spin-%d.png" % (self.magnetic_status, i+1))
            plt.close()

        if spin_n > 1:
            # all spin in on figure
            for i in range(spin_n):
                for j in range(band_min, band_max, 1):
                    if self.efermi == None:
                        plt.plot(self.xcoord_k, [self.eigenval["spin_%d" % (i+1)][k]["energy"][j] for k in range(len(self.eigenval["spin_%d" % (i+1)]))], linewidth=1, label="spin-%d" % (i+1))
                    else:
                        plt.plot(self.xcoord_k, [self.eigenval["spin_%d" % (i+1)][k]["energy"][j] - self.efermi for k in range(len(self.eigenval["spin_%d" % (i+1)]))], linewidth=1, label="spin-%d" % (i+1))
            #plt.title("Band Structure(all spin)")
            plt.xlabel("K")
            plt.ylabel("Energy(eV)")
            plt.xticks(self.locs, self.labels_for_matplotlib)
            plt.grid(b=True, which='major')
            if xrange != None:
                plt.xlim(xmin=xrange[0], xmax=xrange[1])
            if yrange != None:
                plt.ylim(ymin=yrange[0], ymax=yrange[1])            
            plt.savefig("band-structure-%s-spin-all.png" % self.magnetic_status)
            plt.close()     

    def _plot_band_gnuplot(self, bandrange=[0, 1.0], xrange=None, yrange=None):
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

        if self.magnetic_status == "spin-unpolarized":
            spin_n = 1
        elif self.magnetic_status == "spin-polarized":
            spin_n = 2

        nband = len(self.eigenval["spin_1"][0]["energy"])
        band_min = int(bandrange[0] * nband)
        band_max = int(bandrange[1] * nband)
        for i in range(spin_n):
            with open("all-bands-%s-spin-%d.data" % (self.magnetic_status, i+1), 'w') as fout:
                fout.write("# band structure extracted from bands.bs\n")
                if self.efermi == None:
                    fout.write("# efermi: smearing is not used so fermi energy is not calculated in cp2k\n")
                else:
                    fout.write("# efermi: %f (eV) and eigenvalues are in unit of eV\n" % self.efermi)
                for k in range(nband):
                    for j in range(len(self.xcoord_k)):
                        fout.write("%f %f\n" % (self.xcoord_k[j], self.eigenval["spin_%d" % (i+1)][j]["energy"][k]))
                    fout.write("\n")
            with open("specified-bands-%s-spin-%d.data" % (self.magnetic_status, i+1), 'w') as fout:
                fout.write("# band structure extracted from bands.bs\n")
                if self.efermi == None:
                    fout.write("# efermi: smearing is not used so fermi energy is not calculated in cp2k\n")
                else:
                    fout.write("# efermi: %f (eV) and eigenvalues are in unit of eV\n" % self.efermi)
                for k in range(band_min, band_max, 1):
                    for j in range(len(self.xcoord_k)):
                        fout.write("%f %f\n" % (self.xcoord_k[j], self.eigenval["spin_%d" % (i+1)][j]["energy"][k]))
                    fout.write("\n")

            with open("all-bands-%s-spin-%d.gnuplot" % (self.magnetic_status, i+1), 'w') as fout:
                fout.write("set terminal gif\n")
                fout.write("set output 'all-bands-%s-spin-%d.gif'\n" % (self.magnetic_status, i+1))
                fout.write("unset key\n")
                fout.write("set parametric\n")
                #fout.write("set title 'Band Structure (Spin %d)'\n" % (i+1))
                fout.write("set xlabel 'K'\n")
                fout.write("set ylabel 'Energy(eV)'\n")
                fout.write("set xtics(")
                for j in range(len(self.labels_for_gnuplot)-1):
                    fout.write("'%s' %f, " % (self.labels_for_gnuplot[j], self.locs[j]))
                fout.write("'%s' %f)\n" % (self.labels_for_gnuplot[-1], self.locs[-1]))
                fout.write("set grid xtics ytics\n")
                fout.write("set autoscale\n")
                if xrange != None:
                    fout.write("set xrange [%f:%f]\n" % (xrange[0], xrange[1]))
                if yrange != None:
                    fout.write("set yrange [%f:%f]\n" % (yrange[0], yrange[1]))
                if self.efermi == None:
                    fout.write("plot 'all-bands-%s-spin-%d.data' using 1:2 w l\n" % (self.magnetic_status, i+1))
                else:
                    fout.write("plot 'all-bands-%s-spin-%d.data' using 1:($2-%f) w l\n" % (self.magnetic_status, i+1, self.efermi))
                    
                
            os.system("gnuplot all-bands-%s-spin-%d.gnuplot" % (self.magnetic_status, i+1))

            with open("specified-bands-%s-spin-%d.gnuplot" % (self.magnetic_status, i+1), 'w') as fout:
                fout.write("set terminal gif\n")
                fout.write("set output 'specified-bands-%s-spin-%d.gif'\n" % (self.magnetic_status, i+1))
                fout.write("unset key\n")
                fout.write("set parametric\n")
                #fout.write("set title 'Band Structure (Spin %d)'\n" % (i+1))
                fout.write("set xlabel 'K'\n")
                fout.write("set ylabel 'Energy(eV)'\n")
                fout.write("set xtics(")
                for j in range(len(self.labels_for_gnuplot)-1):
                    fout.write("'%s' %f, " % (self.labels_for_gnuplot[j], self.locs[j]))
                fout.write("'%s' %f)\n" % (self.labels_for_gnuplot[-1], self.locs[-1]))
                fout.write("set grid xtics ytics\n")
                fout.write("set autoscale\n")
                if xrange != None:
                    fout.write("set xrange [%f:%f]\n" % (xrange[0], xrange[1]))
                if yrange != None:
                    fout.write("set yrange [%f:%f]\n" % (yrange[0], yrange[1]))                
                if self.efermi == None:
                    fout.write("plot 'specified-bands-%s-spin-%d.data' using 1:2 w l\n" % (self.magnetic_status, i+1))
                else:
                    fout.write("plot 'specified-bands-%s-spin-%d.data' using 1:($2-%f) w l\n" % (self.magnetic_status, i+1, self.efermi))
                    
            os.system("gnuplot specified-bands-%s-spin-%d.gnuplot" % (self.magnetic_status, i+1))
            
    def plot_band(self, bandrange=[0, 1.0], engine="matplotlib", xrange=None, yrange=None):
        """
        :parama engine:
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
        if engine == "matplotlib":
            self._plot_band_matplotlib(bandrange=bandrange, xrange=xrange, yrange=yrange)
        elif engine == "gnuplot":
            self._plot_band_gnuplot(bandrange=bandrange, xrange=xrange, yrange=yrange)
        #

    def print_gap(self):
        """
        print out the gap of the band
        """
        if self.magnetic_status == "spin-unpolarized":
            spin_n = 1
        elif self.magnetic_status == "spin-polarized":
            spin_n = 2
        
        print("===============================================================\n")
        print("            Info about the gap of the system\n")
        print("---------------------------------------------------------------\n")
        if self.efermi == None:
            print("Smearing is not used, so Efermi is not calculated\n")
            print("failed to calc the energy gap\n")
            return
        
        nband = len(self.eigenval["spin_1"][0]["energy"])
        nkpoints = len(self.eigenval["spin_1"])
        if spin_n == 1:
            print("The calculation is done for only one spin\n")
            is_metallic = False
            for i in range(nband):
                band_energies = [self.eigenval["spin_1"][k]["energy"][i] for k in range(nkpoints)]
                if min(band_energies) < self.efermi < max(band_energies):
                    is_metallic = True
                    break
            if is_metallic == True:
                print("The system is metallic!\n")
                return
            # is_metallic == False:
            homo = 0
            lumo = 0
            valence = []
            conduction = []
            for i in range(nband):
                band_energies = [self.eigenval["spin_1"][k]["energy"][i] for k in range(nkpoints)]
                if max(band_energies) <= self.efermi:
                    valence.append(i)
                    continue
                elif min(band_energies) >= self.efermi:
                    conduction.append(i)
                    continue
                else:
                    print("WARNING: band %d can not be classified to valence or conduction\n" % (i))
            homo = max(valence)
            lumo = min(conduction)
            gap = min([self.eigenval["spin_1"][k]["energy"][lumo] for k in range(nkpoints)]) - max([self.eigenval["spin_1"][k]["energy"][homo] for k in range(nkpoints)])
            print("The system is semiconducting or insulating\n")
            print("And the gap is %f eV\n" % gap)
        elif spin_n == 2:
            print("The calculation is done for two spin\n")
            print("Spin 1:\n")
            is_metallic_spin_1 = False
            for i in range(nband):
                band_energies = [self.eigenval["spin_1"][k]["energy"][i] for k in range(nkpoints)]
                if min(band_energies) < self.efermi < max(band_energies):
                    is_metallic_spin_1 = True
                    break

            if is_metallic_spin_1 == True:
                print("The system is metallic in spin 1!\n")
            else:
                # is_metallic == False:
                homo = 0
                lumo = 0
                valence = []
                conduction = []
                for i in range(nband):
                    band_energies = [self.eigenval["spin_1"][k]["energy"][i] for k in range(nkpoints)]
                    if max(band_energies) <= self.efermi:
                        valence.append(i)
                        continue
                    elif min(band_energies) >= self.efermi:
                        conduction.append(i)
                        continue
                    else:
                        print("WARNING: band %d in spin 1 can not be classified to valence or conduction\n" % (i))
                homo = max(valence)
                lumo = min(conduction)
                gap = min([self.eigenval["spin_1"][k]["energy"][lumo] for k in range(nkpoints)]) - max([self.eigenval["spin_1"][k]["energy"][homo] for k in range(nkpoints)])
                print("The system is semiconducting or insulating in spin 1\n")
                print("And the gap in spin 1 is %f eV\n" % gap)
            
            print("Spin 2:\n")
            is_metallic_spin_2 = False
            for i in range(nband):
                band_energies = [self.eigenval["spin_2"][k]["energy"][i] for k in range(nkpoints)]
                if min(band_energies) < self.efermi < max(band_energies):
                    is_metallic_spin_2 = True
                    break

            if is_metallic_spin_2 == True:
                print("The system is metallic in spin 2!\n")
            else:
                # is_metallic == False:
                homo = 0
                lumo = 0
                valence = []
                conduction = []
                for i in range(nband):
                    band_energies = [self.eigenval["spin_2"][k]["energy"][i] for k in range(nkpoints)]
                    if max(band_energies) <= self.efermi:
                        valence.append(i)
                        continue
                    elif min(band_energies) >= self.efermi:
                        conduction.append(i)
                        continue
                    else:
                        print("WARNING: band %d in spin 2 can not be classified to valence or conduction\n" % (i))
                homo = max(valence)
                lumo = min(conduction)
                gap = min([self.eigenval["spin_2"][k]["energy"][lumo] for k in range(nkpoints)]) - max([self.eigenval["spin_2"][k]["energy"][homo] for k in range(nkpoints)])
                print("The system is semiconducting or insulating in spin 2\n")
                print("And the gap in spin 2 is %f eV\n" % gap)


    def print_effective_mass(self):
        """
        print out the effective mass of the band(top of homo and buttom of lumo)
        """
        if self.magnetic_status == "spin-unpolarized":
            spin_n = 1
        elif self.magnetic_status == "spin-polarized":
            spin_n = 2
        
        print("===============================================================\n")
        print("            Info about the effective Mass of the system\n")
        print("---------------------------------------------------------------\n")
        if self.efermi == None:
            print("Smearing is not used, so Efermi is not calculated\n")
            print("failed to calc the effective mass\n")
            return
            
        nband = len(self.eigenval["spin_1"][0]["energy"])
        nkpoints = len(self.eigenval["spin_1"])
        if spin_n == 1:
            print("The calculation is done for only one spin\n")
            is_metallic = False
            for i in range(nband):
                band_energies = [self.eigenval["spin_1"][k]["energy"][i] for k in range(nkpoints)]
                if min(band_energies) < self.efermi < max(band_energies):
                    is_metallic = True
                    break
            if is_metallic == True:
                print("The system is metallic!\n")
                return
            # is_metallic == False:
            homo = 0
            lumo = 0
            valence = []
            conduction = []
            for i in range(nband):
                band_energies = [self.eigenval["spin_1"][k]["energy"][i] for k in range(nkpoints)]
                if max(band_energies) <= self.efermi:
                    valence.append(i)
                    continue
                elif min(band_energies) >= self.efermi:
                    conduction.append(i)
                    continue
                else:
                    print("WARNING: band %d can not be classified to valence or conduction\n" % (i))
            homo = max(valence)
            lumo = min(conduction)
            
            top_energy = max([self.eigenval["spin_1"][k]["energy"][homo] for k in range(nkpoints)])
            low_energy = min([self.eigenval["spin_1"][k]["energy"][lumo] for k in range(nkpoints)])
            for k in range(nkpoints):
                if self.eigenval["spin_1"][k]["energy"][homo] == top_energy:
                    top_peak_k = k
            for k in range(nkpoints):
                if self.eigenval["spin_1"][k]["energy"][lumo] == low_energy:
                    low_valley_k = k
            # now interpolate datas around top_peak_k and low_valley_k
            # on the left we choose num_k_left, on the right we choose num_k_right
            num_k_left = 5 
            num_k_right = 5
            homo_x = [self.xcoord_k[top_peak_k+i] for i in range(-num_k_left, +num_k_right+1)]
            homo_y = [self.eigenval["spin_1"][top_peak_k+i]["energy"][homo] for i in range(-num_k_left, +num_k_right+1)]
            lumo_x = [self.xcoord_k[low_valley_k+i] for i in range(-num_k_left, +num_k_right+1)]
            lumo_y = [self.eigenval["spin_1"][low_valley_k+i]["energy"][lumo] for i in range(-num_k_left, +num_k_right+1)]        
            # unit of y is eV
            # unit of x is 1/a namely 1/Angstrom
            # we now transfor m unit of y to Hartree and unit of x to bohr-1
            Bohr = 0.529177208
            for i in range(len(homo_x)):
                homo_x[i] =  homo_x[i] * Bohr                
            for i in range(len(lumo_x)):
                lumo_x[i] = lumo_x[i] * Bohr
            Hartree = 27.211396
            for i in range(len(homo_y)):
                homo_y[i] = homo_y[i] / Hartree
            for i in range(len(lumo_y)):
                lumo_y[i] = lumo_y[i] / Hartree
                
            homo_result = np.polyfit(homo_x, homo_y, deg=2) # fit to a quadratic function: y = a * x**2 + b * x + c
            lumo_result = np.polyfit(lumo_x, lumo_y, deg=2)
            m_effective_homo = 1 / (2*homo_result[0]) # 1/(2a)
            m_effective_lumo = 1 / (2*lumo_result[0])
            
            print("The system is semiconducting or insulating\n")
            print("And the effective mass at HOMO is %f\n" % (m_effective_homo))
            print("the corresponding E-k is %f %f\n" % (top_energy, self.xcoord_k[top_peak_k]))
            print("And the effective mass at LUMO is %f\n" % (m_effective_lumo))
            print("the corresponding E-k is %f %f\n" % (low_energy, self.xcoord_k[low_valley_k]))
            print("The Efermi is: %f\n" % self.efermi)
            print("The value is not right as the unit of k is not rightly handled now\n")
        elif spin_n == 2:
            print("The calculation is done for two spin\n")
            print("Spin 1:\n")
            is_metallic_spin_1 = False
            for i in range(nband):
                band_energies = [self.eigenval["spin_1"][k]["energy"][i] for k in range(nkpoints)]
                if min(band_energies) < self.efermi < max(band_energies):
                    is_metallic_spin_1 = True
                    break

            if is_metallic_spin_1 == True:
                print("The system is metallic in spin 1!\n")
            else:
                # is_metallic == False:
                homo = 0
                lumo = 0
                valence = []
                conduction = []
                for i in range(nband):
                    band_energies = [self.eigenval["spin_1"][k]["energy"][i] for k in range(nkpoints)]
                    if max(band_energies) <= self.efermi:
                        valence.append(i)
                        continue
                    elif min(band_energies) >= self.efermi:
                        conduction.append(i)
                        continue
                    else:
                        print("WARNING: band %d in spin 1 can not be classified to valence or conduction\n" % (i))
                homo = max(valence)
                lumo = min(conduction)
                
                top_energy = max([self.eigenval["spin_1"][k]["energy"][homo] for k in range(nkpoints)])
                low_energy = min([self.eigenval["spin_1"][k]["energy"][lumo] for k in range(nkpoints)])
                for k in range(nkpoints):
                    if self.eigenval["spin_1"][k]["energy"][homo] == top_energy:
                        top_peak_k = k
                for k in range(nkpoints):
                    if self.eigenval["spin_1"][k]["energy"][lumo] == low_energy:
                        low_valley_k = k
                # now interpolate datas around top_peak_k and low_valley_k
                # on the left we choose num_k_left, on the right we choose num_k_right
                num_k_left = 5 
                num_k_right = 5
                homo_x = [self.xcoord_k[top_peak_k+i] for i in range(-num_k_left, +num_k_right+1)]
                homo_y = [self.eigenval["spin_1"][top_peak_k+i]["energy"][homo] for i in range(-num_k_left, +num_k_right+1)]
                lumo_x = [self.xcoord_k[low_valley_k+i] for i in range(-num_k_left, +num_k_right+1)]
                lumo_y = [self.eigenval["spin_1"][low_valley_k+i]["energy"][lumo] for i in range(-num_k_left, +num_k_right+1)]        
                # unit of y is eV
                # unit of x is 1/a namely 1/Angstrom
                # we now transfor m unit of y to Hartree and unit of x to bohr-1
                Bohr = 0.529177208
                for i in range(len(homo_x)):
                    homo_x[i] =  homo_x[i] * Bohr                
                for i in range(len(lumo_x)):
                    lumo_x[i] = lumo_x[i] * Bohr
                Hartree = 27.211396
                for i in range(len(homo_y)):
                    homo_y[i] = homo_y[i] / Hartree
                for i in range(len(lumo_y)):
                    lumo_y[i] = lumo_y[i] / Hartree
                
                homo_result = np.polyfit(homo_x, homo_y, deg=2) # fit to a quadratic function: y = a * x**2 + b * x + c
                lumo_result = np.polyfit(lumo_x, lumo_y, deg=2)
                m_effective_homo = 1 / (2*homo_result[0]) # 1/(2a)
                m_effective_lumo = 1 / (2*lumo_result[0])
            
                print("The system is semiconducting or insulating in spin 1\n")
                print("And the effective mass at HOMO of spin 1 is %f\n" % (m_effective_homo))
                print("the corresponding E-k is %f %f\n" % (top_energy, self.xcoord_k[top_peak_k]))
                print("And the effective mass at LUMO of spin 1 is %f\n" % (m_effective_lumo))
                print("the corresponding E-k is %f %f\n" % (low_energy, self.xcoord_k[low_valley_k]))
                print("The Efermi is: %f\n" % self.efermi)
                print("The value is not right as the unit of k is not rightly handled now\n")
            print("Spin 2:\n")
            is_metallic_spin_2 = False
            for i in range(nband):
                band_energies = [self.eigenval["spin_2"][k]["energy"][i] for k in range(nkpoints)]
                if min(band_energies) < self.efermi < max(band_energies):
                    is_metallic_spin_2 = True
                    break

            if is_metallic_spin_2 == True:
                print("The system is metallic in spin 2!\n")
            else:
                # is_metallic == False:
                homo = 0
                lumo = 0
                valence = []
                conduction = []
                for i in range(nband):
                    band_energies = [self.eigenval["spin_2"][k]["energy"][i] for k in range(nkpoints)]
                    if max(band_energies) <= self.efermi:
                        valence.append(i)
                        continue
                    elif min(band_energies) >= self.efermi:
                        conduction.append(i)
                        continue
                    else:
                        print("WARNING: band %d in spin 2 can not be classified to valence or conduction\n" % (i))
                homo = max(valence)
                lumo = min(conduction)
                top_energy = max([self.eigenval["spin_1"][k]["energy"][homo] for k in range(nkpoints)])
                low_energy = min([self.eigenval["spin_1"][k]["energy"][lumo] for k in range(nkpoints)])
                for k in range(nkpoints):
                    if self.eigenval["spin_1"][k]["energy"][homo] == top_energy:
                        top_peak_k = k
                for k in range(nkpoints):
                    if self.eigenval["spin_1"][k]["energy"][lumo] == low_energy:
                        low_valley_k = k
                # now interpolate datas around top_peak_k and low_valley_k
                # on the left we choose num_k_left, on the right we choose num_k_right
                num_k_left = 5 
                num_k_right = 5
                homo_x = [self.xcoord_k[top_peak_k+i] for i in range(-num_k_left, +num_k_right+1)]
                homo_y = [self.eigenval["spin_1"][top_peak_k+i]["energy"][homo] for i in range(-num_k_left, +num_k_right+1)]
                lumo_x = [self.xcoord_k[low_valley_k+i] for i in range(-num_k_left, +num_k_right+1)]
                lumo_y = [self.eigenval["spin_1"][low_valley_k+i]["energy"][lumo] for i in range(-num_k_left, +num_k_right+1)]        
                # unit of y is eV
                # unit of x is 1/a namely 1/Angstrom
                # we now transfor m unit of y to Hartree and unit of x to bohr-1
                Bohr = 0.529177208
                for i in range(len(homo_x)):
                    homo_x[i] =  homo_x[i] * Bohr                
                for i in range(len(lumo_x)):
                    lumo_x[i] = lumo_x[i] * Bohr
                Hartree = 27.211396
                for i in range(len(homo_y)):
                    homo_y[i] = homo_y[i] / Hartree
                for i in range(len(lumo_y)):
                    lumo_y[i] = lumo_y[i] / Hartree
                
                homo_result = np.polyfit(homo_x, homo_y, deg=2) # fit to a quadratic function: y = a * x**2 + b * x + c
                lumo_result = np.polyfit(lumo_x, lumo_y, deg=2)
                m_effective_homo = 1 / (2*homo_result[0]) # 1/(2a)
                m_effective_lumo = 1 / (2*lumo_result[0])
            
                print("The system is semiconducting or insulating in spin 2\n")
                print("And the effective mass of spin 2 at HOMO is %f\n" % (m_effective_homo))
                print("the corresponding E-k is %f %f\n" % (top_energy, self.xcoord_k[top_peak_k]))
                print("And the effective mass of spin 2 at LUMO is %f\n" % (m_effective_lumo))
                print("the corresponding E-k is %f %f\n" % (low_energy, self.xcoord_k[low_valley_k]))
                print("The Efermi is: %f\n" % self.efermi)
                print("The value is not right as the unit of k is not rightly handled now\n")

    def export(self, directory="tmp-cp2k-static", bandrange=[0, 1], engine="matplotlib", xrange=None, yrange=None):
        """
        :parama engine:
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
        self.plot_band(engine=engine,  bandrange=bandrange, xrange=xrange, yrange=yrange)
        os.chdir("../../")


class bands_post_old:
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
                fout.write("unset key\n")
                fout.write("set parametric\n")

                fout.write("set title 'Bandstructure'\n")
                #fout.write("set xlabel 'Kpoint number'\n") # equal to nmber of kpoints
                fout.write("set ylabel 'Energy'\n")
                fout.write("unset key\n")
                fout.write("set xtics(")
                for point in specialk:
                    if point["label"] == "GAMMA":
                        fout.write("'%s' %d, " % ("{/symbol G}", point["k-number"]-1))
                    else:
                        fout.write("'%s' %d, " % (point["label"], point["k-number"]-1)) # minus 1, because in gnuplot x start with 0
                fout.write(")\n")
                fout.write("plot for [i=4:%d] '%s.set-1.csv' u 0:i w l \n" % (n_bands+3, bandsfile))
            os.system("gnuplot bandplot.gp")

        elif option == "matplotlib":
            pass
