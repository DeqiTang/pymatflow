#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import numpy as np
import argparse

from pymatflow.cmd.matflow import get_kpath

"""
usage:
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath for phonon band calc")

    parser.add_argument("--kpath-file", type=str,
            help="file to read the kpath for phonon band calc")

    parser.add_argument("-o", "--qpoints-file", type=str, default="matdyn-qpoints.txt",
            help="output file containings qpoints for matdyn")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

    kpath = get_kpath(kpath_manual=args.kpath_manual, kpath_file=args.kpath_file)

    qpoints = []

    qpoints.append([kpath[0][0], kpath[0][1], kpath[0][2], 0.0000000,  kpath[0][3]])
    for i in range(1, len(kpath)):
        if kpath[i-1][4] == "|":
            qpoints.append([kpath[i][0], kpath[i][1], kpath[i][2], qpoints[-1][3], kpath[i][3]])
        else:
            step_kx = (kpath[i][0] - kpath[i-1][0]) / kpath[i-1][4]
            step_ky = (kpath[i][1] - kpath[i-1][1]) / kpath[i-1][4]
            step_kz = (kpath[i][2] - kpath[i-1][2]) / kpath[i-1][4]
            step_xcoord_k = np.sqrt((kpath[i][0]-kpath[i-1][0])**2+(kpath[i][1]-kpath[i-1][1])**2+(kpath[i][2]-kpath[i-1][2])**2) / kpath[i-1][4]
            for j in range(kpath[i-1][4]-1):
                qpoints.append([
                    kpath[i-1][0] + step_kx * (j+1),
                    kpath[i-1][1] + step_ky * (j+1),
                    kpath[i-1][2] + step_kz * (j+1),
                    qpoints[-1][3] + step_xcoord_k,
                    None
                    ])
            qpoints.append([kpath[i][0], kpath[i][1], kpath[i][2], qpoints[-1][3]+step_xcoord_k, kpath[i][3]])


    with open(args.qpoints_file, 'w') as fout:
        fout.write("%d\n" % len(qpoints))
        for qpoint in qpoints:
            if qpoint[4] == None:
                fout.write("%f %f %f %f\n" % (qpoint[0], qpoint[1], qpoint[2], qpoint[3]))
            else:
                fout.write("%f %f %f %f #%s\n" % (qpoint[0], qpoint[1], qpoint[2], qpoint[3], qpoint[4]))
