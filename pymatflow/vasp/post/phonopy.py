
import os


from pymatflow.base.xyz import base_xyz


class phonopy_post:
    def __init__(self):
        self.supercell_n = [1, 1, 1]
        self.mp = [8, 8, 8]


    def get_kpath(self, kpath):
        self.kpath = kpath

    def get_xyz(self, filepath):
        self.xyz = base_xyz()
        self.xyz.get_xyz(filepath)

    def export(self, directory, engine="matplotlib", tmax=1500):
        os.chdir(directory)
        os.system("ls | grep 'POSCAR-' > pos.data")
        disps = []
        with open("pos.data", 'r') as fin:
            for line in fin:
                disps.append(line.split("\n")[0].split("-")[1])
        os.chdir("../")

        os.chdir(directory)
        os.system("mkdir -p post-processing")
        # generate the result analysis bash script and necessary config files

        with open("post-processing/mesh.conf", 'w') as fout:
            fout.write("ATOM_NAME =")
            for element in self.xyz.specie_labels:
                fout.write(" %s" % element)
            fout.write("\n")
            fout.write("DIM = %d %d %d\n" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))
            fout.write("MP = %d %d %d\n" % (self.mp[0], self.mp[1], self.mp[2]))


        with open("post-processing/pdos.conf", 'w') as fout:
            fout.write("ATOM_NAME =")
            for element in self.xyz.specie_labels:
                fout.write(" %s" % element)
            fout.write("\n")
            fout.write("DIM = %d %d %d\n" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))
            fout.write("MP = %d %d %d\n" % (self.mp[0], self.mp[1], self.mp[2]))
            fout.write("PDOS = 1 2, 3 4 5 5\n")

        with open("post-processing/band.conf", 'w') as fout:
            fout.write("ATOM_NAME =")
            for element in self.xyz.specie_labels:
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
            fout.write("DIM = %d %d %d\n" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))
            fout.write("BAND =")
            for i in range(len(self.kpath) - 1):
                if self.kpath[i][4] != "|":
                    fout.write(" %f %f %f" % (self.kpath[i][0], self.kpath[i][1], self.kpath[i][2]))
                else:
                    fout.write(" %f %f %f," % (self.kpath[i][0], self.kpath[i][1], self.kpath[i][2]))
            fout.write(" %f %f %f" % (self.kpath[-1][0], self.kpath[-1][1], self.kpath[-1][2]))
            fout.write("\n")
            fout.write("BAND_LABELS =")
            for i in range(len(self.kpath) - 1):
                if self.kpath[i][4] != "|":
                    if self.kpath[i][3].upper() == "GAMMA":
                        fout.write(" $\Gamma$")
                    else:
                        fout.write(" $%s$" % self.kpath[i][3])
                else:
                    if self.kpath[i][3].upper() == "GAMMA":
                        fout.write(" $\Gamma$,")
                    else:
                        fout.write(" $%s$," % self.kpath[i][3])
            if self.kpath[-1][3].upper() == "GAMMA":
                fout.write(" $\Gamma$")
            else:
                fout.write(" $%s$" % self.kpath[-1][3])
            fout.write("\n")

        with open("post-processing/phonopy-analysis.sh", 'w') as fout:
            fout.write("#!/bin/bash\n\n")
            fout.write("cp ../phonopy_disp.yaml ./\n")
            fout.write("cp ../disp.yaml ./\n")
            # From phonopy v2.0, displacements are written into "phonopy_disp.yaml".
            # "disp.yaml" is still supported for reading except for Wien2k interface, 
            # so we now cp both possible phonopy_disp.yaml and disp.yaml to current dir
            fout.write("cp ../POSCAR ./\n")
            fout.write("# generate the FORCE_SET\n")
            fout.write("phonopy -f ../disp-{001..%s}/vasprun.xml\n" % (disps[-1]))
            fout.write("# plot The density of states (DOS)\n")
            fout.write("phonopy -p mesh.conf -s\n")
            fout.write("# Thermal properties are calculated with the sampling mesh by:\n")
            fout.write("phonopy -t mesh.conf --tmax %f\n" % tmax)
            fout.write("# Thermal properties can be plotted by:\n")
            fout.write("phonopy -t -p mesh.conf -s --tmax %f\n" % tmax)
            fout.write("# calculate Projected DOS and plot it\n")
            fout.write("phonopy -p pdos.conf -s\n")
            fout.write("# plot band structure\n")
            fout.write("phonopy -p band.conf -s\n")

        os.system("cd post-processing; bash phonopy-analysis.sh; cd ../")

        import yaml
        with open("post-processing/band.yaml", 'r') as fin:
            band_yaml = yaml.safe_load(fin)
        
        npath = band_yaml["npath"]
        segment_nqpoint = band_yaml["segment_nqpoint"]
        labels = band_yaml["labels"]
        
        with open("post-processing/band.data", 'w') as fout:
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

