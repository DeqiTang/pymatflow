
import os

from pymatflow.base.xyz import base_xyz


class phonon_post:
    def __init__(self):
        self.mp = [8, 8, 8] # default value
        self.supsercell_n = [1, 1, 1] # default value


    def get_kpath(self, kpath):
        self.kpath = kpath

    def get_xyz(filepath):
        self.xyz = base_xyz()
        self.xyz.get_xyz(filepath)


    def export(self, directory):
        # get the disps information
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
            for qpoint in self.kpath:
                if qpoint[4] == None:
                    fout.write(" %f %f %f" % (qpoint[0], qpoint[1], qpoint[2]))
                elif qpoint[4] == "|":
                    fout.write(" %f %f %f," % (qpoint[0], qpoint[1], qpoint[2]))
                else:
                    pass
            fout.write("\n")
            fout.write("BAND_LABELS =")
            for qpoint in self.kpath:
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

        with open("post-processing/phonon-analysis.sh", 'w') as fout:
            fout.write("#!/bin/bash\n\n")
            #fout.write("cp ../phonopy_disp.yaml ./")
            #fout.write("cp ../POSCAR ./\n")
            fout.write("cp ../POSCAR-original-unitcell ./\n")
            fout.write("# generate the FORCE_SET\n")
            fout.write("phonopy --fc ../vasprun.xml\n")
            #fout.write("phonopy -f ../disp-{001..%s}/vasprun.xml\n" % (disps[-1]))
            fout.write("# plot The density of states (DOS)\n")
            fout.write("phonopy -p mesh.conf -c POSCAR-original-unitcell --readfc\n")
            fout.write("# Thermal properties are calculated with the sampling mesh by:\n")
            fout.write("phonopy -t mesh.conf -c POSCAR-original-unitcell --readfc\n")
            fout.write("# Thermal properties can be plotted by:\n")
            fout.write("phonopy -t -p mesh.conf -c POSCAR-original-unitcell --readfc\n")
            fout.write("# calculate Projected DOS and plot it\n")
            fout.write("phonopy -p pdos.conf -c POSCAR-original-unitcell --readfc\n")
            fout.write("# plot band structure\n")
            fout.write("phonopy -p band.conf -c POSCAR-original-unitcell --readfc\n")

        os.system("cd post-processing; bash phonon-analysis.sh; cd ../")

        os.chdir("../")
