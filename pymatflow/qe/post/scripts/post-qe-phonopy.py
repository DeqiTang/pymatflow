#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.qe.base.arts import qe_arts

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of phonopy running", type=str, default="tmp-qe-phonopy")
    parser.add_argument("-f", "--file", help="input structure file", type=str, default=None)
   

    parser.add_argument("--qpoints", type=str, nargs="+", default=None,
            help="manual input qpoints with labels to set the BAND for phonopy analysis")
    
    parser.add_argument("--qpoints-file", type=str, default="kpath-from-seekpath.txt",
            help="file to read the qpoints to set the BAND for phonopy analysis")

    parser.add_argument("--mp", type=int, nargs="+",
            default=[8, 8, 8],
            help="MP for mesh.conf and pdos.conf")

    parser.add_argument("--supercell-n", type=int, nargs="+",
            default=None,
            help="supercell_n used in phonopy calculation")

    args = parser.parse_args()
    
    arts = qe_arts()
    arts.xyz.get_xyz(args.file)
  
    # obtain the qpoints 
    qpoints = [] # [[kx, ky, kz, label, end_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', None], ...]
    # if end_indicator in a qpoint is None, then it will connect to the following point,
    # if end_indicator in a qpoint is '|', then it will not connect to the following point,
    if args.qpoints != None:
        # qpoints from script argument args.qpoints
        for qpoint in args.qpoints:
            if len(qpoint.split()) == 4:
                qpoints.append([
                    float(qpoint.split()[0]),
                    float(qpoint.split()[1]),
                    float(qpoint.split()[2]),
                    qpoint.split()[3].upper(),
                    None,
                    ])
            elif len(qpoint.split()) == 5 and qpoint.split()[4] == "|":
                qpoints.append([
                    float(qpoint.split()[0]),
                    float(qpoint.split()[1]),
                    float(qpoint.split()[2]),
                    qpoint.split()[3].upper(),
                    "|",
                    ])
    else:
        # qpoints read from file specified by args.qpoints_file
        # file is in format like this
        """
        3
        0.0 0.0 0.0 #GAMMA
        x.x x.x x.x #XXX |
        x.x x.x x.x #XXX
        """
        # if there is a '|' behind the label it means the path is 
        # broken after that point!!!
        with open(args.qpoints_file, 'r') as fin:
            qpoints_file = fin.readlines()
        nq = int(qpoints_file[0])
        for i in range(nq):
            if len(qpoints_file[i+1].split()) == 4:
                qpoints.append([
                    float(qpoints_file[i+1].split()[0]),
                    float(qpoints_file[i+1].split()[1]),
                    float(qpoints_file[i+1].split()[2]),
                    qpoints_file[i+1].split()[3].split("#")[1].upper(),
                    None,
                    ])
            elif len(qpoints_file[i+1].split()) == 5 and qpoints_file[i+1].split("\n")[0].split()[4] == "|":
                qpoints.append([
                    float(qpoints_file[i+1].split()[0]),
                    float(qpoints_file[i+1].split()[1]),
                    float(qpoints_file[i+1].split()[2]),
                    qpoints_file[i+1].split()[3].split("#")[1].upper(),
                    '|',
                    ])

    # get the disps information
    os.chdir(args.directory)
    os.system("ls | grep 'supercell-' > pos.data")
    disps = []
    with open("pos.data", 'r') as fin:
        for line in fin:
            disps.append(line.split(".")[0].split("-")[1])
    os.chdir("../") 

    os.chdir(args.directory)

    # generate the result analysis bash script and necessary config files

    with open("mesh.conf", 'w') as fout:
        fout.write("ATOM_NAME =")
        for element in arts.xyz.specie_labels:
            fout.write(" %s" % element)
        fout.write("\n")
        fout.write("DIM = %d %d %d\n" % (args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))
        fout.write("MP = %d %d %d\n" % (args.mp[0], args.mp[1], args.mp[2]))
 
            
                        
    with open("pdos.conf", 'w') as fout:
        fout.write("ATOM_NAME =")
        for element in arts.xyz.specie_labels:
            fout.write(" %s" % element)
        fout.write("\n")
        fout.write("DIM = %d %d %d\n" % (args.supercell_n[0], args.supercell_n[1], args.supercell_n[2]))
        fout.write("MP = %d %d %d\n" % (args.mp[0], args.mp[1], args.mp[2]))
        fout.write("PDOS = 1 2, 3 4 5 5\n")

    with open("band.conf", 'w') as fout:
        fout.write("ATOM_NAME =")
        for element in arts.xyz.specie_labels:
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
        for qpoint in qpoints:
            if qpoint[4] == None:
                fout.write(" %f %f %f" % (qpoint[0], qpoint[1], qpoint[2]))
            elif qpoint[4] == "|":
                fout.write(" %f %f %f," % (qpoint[0], qpoint[1], qpoint[2]))
            else:
                pass
        fout.write("\n")
        fout.write("BAND_LABELS =")
        for qpoint in qpoints:
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

    with open("phonopy-analysis.sh", 'w') as fout:
        fout.write("#!/bin/bash\n\n")
        fout.write("# generate the FORCE_SET\n")
        fout.write("phonopy --qe -f supercell-{001..%s}.out\n" % (disps[-1]))
        fout.write("# plot The density of states (DOS)\n")
        fout.write("phonopy --qe -p mesh.conf -c pos.in\n")
        fout.write("# Thermal properties are calculated with the sampling mesh by:\n")
        fout.write("phonopy --qe -t mesh.conf -c pos.in\n")
        fout.write("# Thermal properties can be plotted by:\n")
        fout.write("phonopy --qe -t -p mesh.conf -c pos.in\n")
        fout.write("# calculate Projected DOS and plot it\n")
        fout.write("phonopy --qe -p pdos.conf -c pos.in\n")
        fout.write("# plot band structure\n")
        fout.write("phonopy --qe -p band.conf -c pos.in\n")

    os.system("bash phonopy-analysis.sh")

    os.chdir("../")

