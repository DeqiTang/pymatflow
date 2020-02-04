#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse


"""
usage:
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-static",
            help="Directory for the static running.")

    parser.add_argument("--bands-out", type=str, default="bands.out",
            help="bands.x output file")

    parser.add_argument("--bands-dat-gnu", type=str, default="bands.dat.gnu",
            help="bands.x output band structure file for gnuplot")

    parser.add_argument("--qpoints-file", type=str, default="matdyn-qpoints.txt",
            help="output file containings qpoints for matdyn")

    parser.add_argument("--pwx-bands-in", type=str, default="static-bands.in",
            help="input for pw.x bands calculation")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

    
    kpoints = []  
    # kpoints = [[kx, ky, kz, xcoord, label], ...] like [[0.0, 0,0, 0.0, 0.0, 'GAMMA'], ...]
    # if label in kpoint in kpoints is None, then kpoint is not a special k
    # else it is a specialk, with label be the special k label.
    # kx, ky, kz are from bands.x output ,it seems to be in tpiba_b type
    # in unit of 2pi/a

    # get all the k points used in calculation.
    with open(os.path.join(args.directory, args.bands_out), 'r') as fin:
        #bands_out = fin.readlines()
        for line in fin:
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "xk=(":
                kpoints.append([
                    float(line.split()[1].split(",")[0]),
                    float(line.split()[2].split(",")[0]),
                    float(line.split()[3]),
                    None, # x coordinates read from bands.dat.gnu file
                    None, # this information is later set by specialk
                    ])

    # get the cooresponding x coordinates for plot for all the k points used in calculation.
    with open(os.path.join(args.directory, args.bands_dat_gnu), 'r') as fin:
        i = 0
        for line in fin:
            if len(line.split()) == 0:
                break
            kpoints[i][3] = float(line.split()[0])
            i = i + 1


    # --------------------------------------------------------------------------
    # get the high symmetry k points label from pw.x bands calculation
    # specialk: [{'label': 'GAMMA', 'coord': [kx, ky, kz], 'xcoord': float}, ...]
    # type of coord in specialk is in type specified in pw.x input file for bands
    # calculation, it can be tpiba_ba or crystal_b(more likely from my point of 
    # view), however, it doesn't matter as they will all by converted to tpiba_b
    # type in the output bands.x calculation(this is not for sure yet!), so the
    # qpoints for matdyn.x are all dependent on the k point from the bands.x
    # output ,and we only use pw.x bands calculation input file to get information
    # about labels of high symmetry k point.
    specialk = []
    with open(os.path.join(args.directory, args.pwx_bands_in), 'r') as fin:
        pwxbandsin = fin.readlines()
    nspecialk = 0
    special_k_begin = 0
    special_k_end = 0
    for i in range(len(pwxbandsin)):
        if len(pwxbandsin[i].split()) == 0:
            continue
        if pwxbandsin[i].split()[0] == "K_POINTS":
            nspecialk = int(pwxbandsin[i+1].split()[0])
            special_k_begin = i + 2
            special_k_end = i + 1 + nspecialk
        
    for i in range(special_k_begin, special_k_end + 1):
        print(pwxbandsin[i].split("#"))
        kpoint = {"label": pwxbandsin[i].split("#")[1].split()[0], "coord": [float(pwxbandsin[i].split()[0]), float(pwxbandsin[i].split()[1]), float(pwxbandsin[i].split()[2])], "xcoord": None}
        specialk.append(kpoint)

    # get the x coordinate of the high symmetry k point from band.x output file
    with open(os.path.join(args.directory, args.bands_out), 'r') as fin:
        bandsxout = fin.readlines()
    xcoord_begin = 0
    for i in range(len(bandsxout)):
        if len(bandsxout[i].split()) == 0:
            continue
        if bandsxout[i].split()[0] == "high-symmetry":
            xcoord_begin = i
            break # break for loop in the first high-symmetry
    for i in range(nspecialk):
        print(bandsxout[xcoord_begin+i])
        specialk[i]["xcoord"] = float(bandsxout[xcoord_begin+i].split()[7])
    #
    # ----------------------------------------------------------------------------

    # transfer information from specialk to kpoints
    for kpoint in kpoints:
        for point in specialk:
            if abs(kpoint[3] - point["xcoord"]) < 1.0e-4: # 1.0e-4 can be an appropriate value to judge the equal of the xcoord
                kpoint[4] = point["label"]
                break
   
    # write the  qpoints file
    with open(args.qpoints_file, 'w') as fout:
        fout.write("%d\n" % len(kpoints))
        for kpoint in kpoints:
            if kpoint[4] == None:
                fout.write("%f %f %f %f\n" % (kpoint[0], kpoint[1], kpoint[2], kpoint[3]))
            else:
                fout.write("%f %f %f %f #%s\n" % (kpoint[0], kpoint[1], kpoint[2], kpoint[3], kpoint[4]))
        #
