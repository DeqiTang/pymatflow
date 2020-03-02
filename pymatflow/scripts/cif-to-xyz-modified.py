#!/usr/bin/env python

import os
import argparse


"""
rely on cif2cell: pip install cif2cell
in order to install cif2cell you must make sure python-dev
or python-devel(yum) is installed
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--cif", type=str,
            help="the cif structure file")
    parser.add_argument("-o", "--xyz", type=str,
            help="the output xyz file")
    args = parser.parse_args()

    # output inofrmation
    print("=====================================\n")
    print("      cif-to-xyz-modified.xyz\n")
    print("-------------------------------------\n")
    print("convert from cif to xyz-modified\n")
    print("with the help from cif2cell\n")

    os.system("cif2cell %s --cartesian > tmp-cif2cell-xxx.dat" % args.cif)
    with open("tmp-cif2cell-xxx.dat", 'r') as fin:
        lines = fin.readlines()
    os.system("rm tmp-cif2cell-xxx.dat")
    # get the begin line number of Bravis lattice vectors:
    for i in range(len(lines)):
        if len(lines[i].split()) == 0:
            continue
        if lines[i].split()[0] == "Bravais":
            lat_vec_begin = i + 1
    # get the begin and end line number of All sites for atoms
    for i in range(len(lines)):
        if len(lines[i].split()) == 0:
            continue
        if lines[i].split()[0] == "All" and lines[i].split()[1] == "sites,":
            atomic_coord_begin = i + 2
    for i in range(atomic_coord_begin, len(lines)):
        if len(lines[i].split()) == 0:
            atomic_coord_end = i - 1
            break
            # must break here
    # generate the modifed version of xyz file
    with open(args.xyz, 'w') as fout:
        fout.write("%d\n" % (atomic_coord_end - atomic_coord_begin + 1))
        # write cell
        #fout.write("cell: %f %f %f | %f %f %f | %f %f %f\n" % (cell[0]0), cell[0][1], cell[0][2], cell[1][0], cell[1][1], cell[1][2], cell[2][0], cell[2][1], cell[2][2])
        fout.write("cell: ")
        fout.write(lines[lat_vec_begin].split("\n")[0])
        fout.write(" | ")
        fout.write(lines[lat_vec_begin+1].split("\n")[0])
        fout.write(" | ")
        fout.write(lines[lat_vec_begin+2].split("\n")[0])
        fout.write("\n")
        # wirte coordinates
        for i in range(atomic_coord_begin, atomic_coord_end + 1):
            fout.write(lines[i])
    #
