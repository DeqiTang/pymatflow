"""
post_bands:
    now post_bands can extract band data from vasprun.xml, but it suppose
    the kpoints in vasprun.xml are in crystal coordinates. and it will build
    the kpoints length: xcoord_k from the high symmetry line and the corresponding
    basis for reciprocal space.

    actually in vasp reci_basis is the reciprocal of the real space basis, namely
    b1 = 1 / a1, b2 = 1 / a2 and b3 = 1 / a3.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from xml.etree.ElementTree import parse

class post_bands:
    def __init__(self):
        pass

    def get_efermi(self, vasprun="vasprun.xml"):
        """
        we set efermi in an individual function because we can choose to get efermi from nscf run
        or scf run in this way.
        if you want to get efermi from the scf run specify the vasprun.xml for the scf
        if you want to get efermi from the nscf run specify the vasprun.xml for the nscf
        """
        vasprun_xml = parse(vasprun)
        self.efermi = float(vasprun_xml.getroot().find("calculation").find("dos").find("i").text)    

    def get_kpath_and_vasprun(self, kpath, vasprun="vasprun.xml"):
        # vasprun.xml
        self.vasprun = parse(vasprun)

        # kpath
        self.kpath = kpath

        self.get_xcoord_k()
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

    def get_xcoord_k(self):
        """
        Note:
            xcoord_k is the x axis of the band structure plot
            let's see how it is build from kpoints and the
            crystal lattice or reciprocal lattice.
        """
        self.xcoord_k = []
        # get the lattice parameter and the reciprocal basis
        a1 = np.sqrt(
                float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[0].getchildren()[0].text.split()[0])**2
                + float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[0].getchildren()[0].text.split()[1])**2
                + float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[0].getchildren()[0].text.split()[2])**2
                )
        a2 = np.sqrt(
                float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[0].getchildren()[1].text.split()[0])**2
                + float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[0].getchildren()[1].text.split()[1])**2
                + float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[0].getchildren()[1].text.split()[2])**2
                )
        a3 = np.sqrt(
                float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[0].getchildren()[2].text.split()[0])**2
                + float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[0].getchildren()[2].text.split()[1])**2
                + float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[0].getchildren()[2].text.split()[2])**2
                )

        b1 = np.sqrt(
                float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[2].getchildren()[0].text.split()[0])**2
                + float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[2].getchildren()[0].text.split()[1])**2
                + float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[2].getchildren()[0].text.split()[2])**2
                )

        b2 = np.sqrt(
                float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[2].getchildren()[1].text.split()[0])**2
                + float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[2].getchildren()[1].text.split()[1])**2
                + float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[2].getchildren()[1].text.split()[2])**2
                )

        b3 = np.sqrt(
                float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[2].getchildren()[2].text.split()[0])**2
                + float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[2].getchildren()[2].text.split()[1])**2
                + float(self.vasprun.findall("structure")[0].getchildren()[0].getchildren()[2].getchildren()[2].text.split()[2])**2
                )

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
                        
        # the old code to get kpoints from vasprun.xml and build xcoord_k
        #divisions = int(self.vasprun.find("kpoints").find("generation").find("i").text)

        #n_segment = int(len(self.vasprun.find("kpoints").findall("varray")[0].getchildren()) / divisions)

        #for i in range(n_segment):
            # xcoord_k of the first point of each segment equals to the
            # last xcoord_k of the previous segment
        #    if i == 0:
        #        # begin from 0.0000000
        #        self.xcoord_k.append(0.00000000)
        #    else:
        #        self.xcoord_k.append(self.xcoord_k[-1])

            # the step in the xcoord_k for each segment is different and it is actually
            # the distance between the two high symmetry kpoint in unit of reciprocal coordinates
            # divided by (divisions - 1)
        #    step = 0
        #    delta_b_1 = b1*float(self.vasprun.find("kpoints").findall("varray")[0].getchildren()[i*divisions].text.split()[0]) - b1*float(self.vasprun.find("kpoints").findall("varray")[0].getchildren()[i*15+divisions-1].text.split()[0])

        #    delta_b_2 = b2*float(self.vasprun.find("kpoints").findall("varray")[0].getchildren()[i*divisions].text.split()[1]) - b2*float(self.vasprun.find("kpoints").findall("varray")[0].getchildren()[i*15+divisions-1].text.split()[1])

        #    delta_b_3 = b3*float(self.vasprun.find("kpoints").findall("varray")[0].getchildren()[i*divisions].text.split()[2]) - b3*float(self.vasprun.find("kpoints").findall("varray")[0].getchildren()[i*15+divisions-1].text.split()[2])
        #    step = np.sqrt(delta_b_1**2+delta_b_2**2+delta_b_3**2) / (divisions-1)

        #    for j in range(divisions-1):
        #        self.xcoord_k.append(self.xcoord_k[-1]+step)

    def get_eigenval(self):
        """

        """
        # first check the magnetic_status
        for i in range(len(self.vasprun.getroot().find("parameters"))):
            if self.vasprun.getroot().find("parameters").getchildren()[i].attrib["name"] == "electronic":
                for item in self.vasprun.getroot().find("parameters")[i]:
                    if item.attrib["name"] == "electronic spin":
                        for spin_item in item:
                            if spin_item.attrib["name"] == "ISPIN":
                                ispin = int(spin_item.text.split()[0])
                            elif spin_item.attrib["name"] == "LNONCOLLINEAR":
                                lnoncollinear = spin_item.text.split()[0]
                            elif spin_item.attrib["name"] == "LSORBIT":
                                lsorbit = spin_item.text.split()[0]
        if lsorbit == "T":
            self.magnetic_status = "soc-ispin-%d" % ispin # soc-ispin-1 or soc-ispin-2
            # actually soc-ispin-2 never exists !!!
            # I found when soc and ISPIN=2 are set in the INCAR at the same time,
            # you can find in <incar> section of vasprun.xml ISPIN=2, but it will be turned to 1
            # in <parameters> -> <separator name="electronic"> - > <separator name="electronic spin"> 
            # So VASP will automatically set ISPIN to 1 even when you set ISPIN to 2 
            # if you are considering soc
            # I mean you can set soc and ISPIN=2 at the same time in INCAR, but 
            # VASP will turn ISPIN to 1, and this post script will read it to be ISPIN = 1
            # so never will therere be soc-ispin-2 even when you set it in INCAR.
            # in vasprun.xml, we read the acutally used ISPIN from <parameters>... rather than
            # the input from <incar>.
        else:
            self.magnetic_status = "non-soc-ispin-%d" % ispin # non-soc-ispin-1 or non-soc-ispin-2


        # get the eigenval
        """
        in case of non-soc-ispin-1 or soc-ispin-1
        self.eigenval = {
            "spin_1": [
                {
                    "energy": [],
                    "occupation": [],
                },
                {
                    "energy": [],
                    "occupation": [],
                }
                ...... 
            ], # length equals to number of kpoints
        }

        in case of non-soc-ispin-2 or soc-ispin-2(this actually does not exist)
        self.eigenval = {
            "spin_1": [
                {
                    "energy": [],
                    "occupation": [],
                },
                {
                    "energy": [],
                    "occupation": [],
                }
                ...... 
            ], # length equals to number of kpoints
            "spin_2": [
                {
                    "energy": [],
                    "occupation": [],
                },
                {
                    "energy": [],
                    "occupation": [],
                }
                ...... 
            ], # length equals to number of kpoints        
        }        
        """
        if self.magnetic_status == "non-soc-ispin-1":
            spin_n = 1
        elif self.magnetic_status == "non-soc-ispin-2":
            spin_n = 2
        elif self.magnetic_status == "soc-ispin-1":
            spin_n = 1
        elif self.magnetic_status == "soc-ispin-2":
            spin_n = 2
            # actually soc-ispin-2 never exists !!!
            # I found when soc and ISPIN=2 are set in the INCAR at the same time,
            # you can find in <incar> section of vasprun.xml ISPIN=2, but it will be turned to 1
            # in <parameters> -> <separator name="electronic"> - > <separator name="electronic spin"> 
            # So VASP will automatically set ISPIN to 1 even when you set ISPIN to 2 
            # if you are considering soc
            # I mean you can set soc and ISPIN=2 at the same time in INCAR, but 
            # VASP will turn ISPIN to 1, and this post script will read it to be ISPIN = 1
            # so never will therere be soc-ispin-2 even when you set it in INCAR.
            # in vasprun.xml, we read the acutally used ISPIN from <parameters>... rather than
            # the input from <incar>.

        self.eigenval = {}
        for i in range(spin_n):
            eigenval = []
            for kpoint in self.vasprun.getroot().find("calculation").find("eigenvalues").find("array").findall("set")[0].findall("set")[i]: # notice the .findall("set")[i] indicate which spin data to read
                eigenval.append({})
                eigenval[-1]["energy"] = []
                eigenval[-1]["occupation"] = []
                for j in range(len(kpoint)):
                    eigenval[-1]["energy"].append(float(kpoint[j].text.split()[0]))
                    eigenval[-1]["occupation"].append(float(kpoint[j].text.split()[1]))
            # if HSE is used, in pymatflow the kpoints in the KPOINTS file include the Monkhorst-pack from SCF and 
            # kpoints from kpath. so self.eigeval contains data for both the two part of kpoints.
            # here we remove the Monkhorst-pack part kpoint from self.eigenval.
            nk_scf = len(eigenval) -  len(self.xcoord_k)
            for k in range(nk_scf):
                eigenval[k] = None
            for k in range(nk_scf):
                eigenval.remove(None)
            self.eigenval["spin_%d" % (i+1)] = eigenval

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

        if self.magnetic_status == "non-soc-ispin-1":
            spin_n = 1
        elif self.magnetic_status == "non-soc-ispin-2":
            spin_n = 2
        elif self.magnetic_status == "soc-ispin-1":
            spin_n = 1
        elif self.magnetic_status == "soc-ispin-2":
            spin_n = 2
            # actually soc-ispin-2 never exists !!!
            # I found when soc and ISPIN=2 are set in the INCAR at the same time,
            # you can find in <incar> section of vasprun.xml ISPIN=2, but it will be turned to 1
            # in <parameters> -> <separator name="electronic"> - > <separator name="electronic spin"> 
            # So VASP will automatically set ISPIN to 1 even when you set ISPIN to 2 
            # if you are considering soc
            # I mean you can set soc and ISPIN=2 at the same time in INCAR, but 
            # VASP will turn ISPIN to 1, and this post script will read it to be ISPIN = 1
            # so never will therere be soc-ispin-2 even when you set it in INCAR.
            # in vasprun.xml, we read the acutally used ISPIN from <parameters>... rather than
            # the input from <incar>.

        nband = len(self.eigenval["spin_1"][0]["energy"])
        band_min = int(bandrange[0] * nband)
        band_max = int(bandrange[1] * nband)

        # export the data in gnuplot format
        for i in range(spin_n):
            with open("all-bands-%s-spin-%d.data" % (self.magnetic_status, i+1), 'w') as fout:
                fout.write("# band structure extracted from vasprun.xml\n")
                fout.write("# efermi: %f\n" % self.efermi)
                for k in range(nband):
                    for j in range(len(self.xcoord_k)):
                        fout.write("%f %f\n" % (self.xcoord_k[j], self.eigenval["spin_%d" % (i+1)][j]["energy"][k]))
                    fout.write("\n")
            with open("specified-bands-%s-spin-%d.data" % (self.magnetic_status, i+1), 'w') as fout:
                fout.write("# band structure extracted from vasprun.xml\n")
                fout.write("# efermi: %f\n" % self.efermi)
                for k in range(band_min, band_max, 1):
                    for j in range(len(self.xcoord_k)):
                        fout.write("%f %f\n" % (self.xcoord_k[j], self.eigenval["spin_%d" % (i+1)][j]["energy"][k]))
                    fout.write("\n")
        # make the plot
        for i in range(spin_n):
            for j in range(band_min, band_max, 1):
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

        if self.magnetic_status == "non-soc-ispin-1":
            spin_n = 1
        elif self.magnetic_status == "non-soc-ispin-2":
            spin_n = 2
        elif self.magnetic_status == "soc-ispin-1":
            spin_n = 1
        elif self.magnetic_status == "soc-ispin-2":
            spin_n = 2
            # actually soc-ispin-2 never exists !!!
            # I found when soc and ISPIN=2 are set in the INCAR at the same time,
            # you can find in <incar> section of vasprun.xml ISPIN=2, but it will be turned to 1
            # in <parameters> -> <separator name="electronic"> - > <separator name="electronic spin"> 
            # So VASP will automatically set ISPIN to 1 even when you set ISPIN to 2 
            # if you are considering soc
            # I mean you can set soc and ISPIN=2 at the same time in INCAR, but 
            # VASP will turn ISPIN to 1, and this post script will read it to be ISPIN = 1
            # so never will therere be soc-ispin-2 even when you set it in INCAR.
            # in vasprun.xml, we read the acutally used ISPIN from <parameters>... rather than
            # the input from <incar>.


        nband = len(self.eigenval["spin_1"][0]["energy"])
        band_min = int(bandrange[0] * nband)
        band_max = int(bandrange[1] * nband)
        for i in range(spin_n):
            with open("all-bands-%s-spin-%d.data" % (self.magnetic_status, i+1), 'w') as fout:
                fout.write("# band structure extracted from vasprun.xml\n")
                fout.write("# efermi: %f\n" % self.efermi)
                for k in range(nband):
                    for j in range(len(self.xcoord_k)):
                        fout.write("%f %f\n" % (self.xcoord_k[j], self.eigenval["spin_%d" % (i+1)][j]["energy"][k]))
                    fout.write("\n")
            with open("specified-bands-%s-spin-%d.data" % (self.magnetic_status, i+1), 'w') as fout:
                fout.write("# band structure extracted from vasprun.xml\n")
                fout.write("# efermi: %f\n" % self.efermi)
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
        if self.magnetic_status == "non-soc-ispin-1":
            spin_n = 1
        elif self.magnetic_status == "non-soc-ispin-2":
            spin_n = 2
        elif self.magnetic_status == "soc-ispin-1":
            spin_n = 1
        elif self.magnetic_status == "soc-ispin-2":
            spin_n = 2
        
        print("===============================================================\n")
        print("            Info about the gap of the system\n")
        print("---------------------------------------------------------------\n")
        
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
            print("And the gap is %f\n" % gap)
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
                print("And the gap in spin 1 is %f\n" % gap)
            
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
                print("And the gap in spin 2 is %f\n" % gap)


    def print_effective_mass(self):
        """
        print out the effective mass of the band(top of homo and buttom of lumo)
        """
        if self.magnetic_status == "non-soc-ispin-1":
            spin_n = 1
        elif self.magnetic_status == "non-soc-ispin-2":
            spin_n = 2
        elif self.magnetic_status == "soc-ispin-1":
            spin_n = 1
        elif self.magnetic_status == "soc-ispin-2":
            spin_n = 2
        
        print("===============================================================\n")
        print("            Info about the effective Mass of the system\n")
        print("---------------------------------------------------------------\n")
        
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
                    
                    
    def export(self, directory="tmp-vasp-static", bandrange=[0, 1], engine="matplotlib", xrange=None, yrange=None):
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
