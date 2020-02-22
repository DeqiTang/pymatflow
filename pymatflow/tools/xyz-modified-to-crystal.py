#!/usr/bin/env python

import numpy as np
from pymatflow.base.xyz import base_xyz
import argparse

"""
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--xyz", type=str,
            help="the xyz structure file(modified version)")
    parser.add_argument("-o", "--crystal", type=str,
            help="the structure file in crystal(fractional) coordinates")
    parser.add_argument("--format", type=str, default='qe',
            choices=['qe'], # currently only support qe format
            help="crystal(fractional) coordinate in format for specific software")
    
    args = parser.parse_args()

    xyz = base_xyz()
    xyz.get_xyz(args.xyz)

    print("==========================================================\n")
    print("convert from carteisian(xyz) to crystal(fractional) coords\n")
    print("==========================================================\n")
    print("-------------------------------------------\n")

    if args.format == "qe":
        with open(args.crystal, 'w') as fout:
            # crystal namely fractional coordinate can be convert from cartesian coordinates
            # the conversion process is like transformation of presentation in quantum mechanics
            # the convmat is bulid to do the conversion
            #latcell = np.array(self.xyz.cell)
            #latcell = latcell.reshape(3, 3)
            latcell = np.array(xyz.cell)
            convmat = np.linalg.inv(latcell.T)
            crystal_coord = np.zeros([xyz.natom, 3])
            for i in range(xyz.natom):
                crystal_coord[i] = convmat.dot(np.array([xyz.atoms[i].x, xyz.atoms[i].y, xyz.atoms[i].z]))
            #
            fout.write("ATOMIC_POSITIONS crystal\n")
            for k in range(xyz.natom):
                fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))
            # end crystal type ATOMIC_POSITIONS

    elif args.format == "siesta":
        pass
    elif args.format == "cp2k":
        pass
    else:
        pass
