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

    args = parser.parse_args()

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

    os.chdir("../")
