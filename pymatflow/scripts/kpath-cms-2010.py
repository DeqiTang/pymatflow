#!/usr/bin/env python

import argparse


"""
Reference:
    High-throughput electronic band structure calculations: Challenges and tools
    by Setyawan, Wahyu and Curtarolo, Stefano
    http://dx.doi.org/10.1016/j.commatsci.2010.05.010
"""

def kpath_from_cms_2010(system):
    """
    :param system: crystal system of the cell
    :return kpath
    Note:
        Triclinic: 
            aP
        Monoclinic: 
            mP mS
        Orthorhombic: 
            oP oS oI oF
        Tetragonal: 
            tP tI
        Hexagonal:
            Rhombohedral: hR
            Hexagonal: hP
        Cubic:
            cP cI cF
    """
    if system == "cP":
        # simple cubic
        kpath = []
        kpath.append([0, 0, 0, "GAMMA", 15])
        kpath.append([0.0, 1/2, 0.0, "X", 15])
        kpath.append([1/2, 1/2, 0.0, "M", 15])
        kpath.append([0.0, 0.0, 0.0, "GAMMA", 15])
        kpath.append([1/2, 1/2, 1/2, "R", 15])
        kpath.append([0.0, 1/2, 0.0, "X", "|"])
        kpath.append([1/2, 1/2, 0.0, "M", 15])
        kpath.append([1/2, 1/2, 1/2, "R", "|"])
        return kpath
    if system == "cF":
        # face center cubic
        kpath = []
        kpath.append([0, 0, 0, "GAMMA", 15])
        kpath.append([1/2, 0, 1/2, "X", 15])
        kpath.append([1/2, 1/4, 3/4, "W", 15])
        kpath.append([3/8, 3/8, 3/4, "K", 15])
        kpath.append([0, 0, 0, "GAMMA", 15])
        kpath.append([1/2, 1/2, 1/2, "L", 15])
        kpath.append([5/8, 1/4, 5/8, "U", 15])
        kpath.append([1/2, 1/4, 3/4, "W", 15])
        kpath.append([1/2, 1/2, 1/2, "L", 15])
        kpath.append([3/8, 3/8, 3/4, "K", "|"])
        kpath.append([5/8, 1/4, 5/8, "U", 15])
        kpath.append([1/2, 0, 1/2, "X", "|"])
        return kpath
    if system == "cI":
        # body center cubic
        kpath = []
        kpath.append([0, 0, 0, "GAMMA", 15])
        kpath.append([1/2, -1/2, 1/2, "H", 15])
        kpath.append([0, 0, 1/2, "N", 15])
        kpath.append([0, 0, 0, "GAMMA", 15])
        kpath.append([1/4, 1/4, 1/4, "P", 15])
        kpath.append([1/2, -1/2, 1/2, "H", "|"])
        kpath.append([1/4, 1/4, 1/4, "P", 15])
        kpath.append([0, 0, 1/2, "N", "|"])
        return kpath
    if system == "hP":
        # Hexagonal
        kpath = []
        return kpath
    if system == "hR":
        # Rhombohedral
        kpath = []
        return kpath
    if system == "tP":
        # simple tetragonal
        kpath = []
        return kpath
    if system == "tI":
        # body center tetragonal
        kpath = []
        return kpath
    if system == "oP":
        # simple orthorhombic
        kpath = []
        return kpath
    if system == "oS":
        # base center orthorhombic
        kpath = []
        return kpath
    if system == "oI":
        # body center orthorhombic
        kpath = []
        return kpath
    if system == "oF":
        # face center orthorhombic
        kpath = []
        return kpath
    if system == "mP":
        # simple monoclinic
        kpath = []
        return kpath
    if system == "mS":
        # base center monoclinic
        kpath = []
        return kpath
    if system == "aP":
        # triclinic
        kpath = []
        return kpath



def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--crystal-system", type=str, default=None,
            choices=["cP", "cI", "cF", "hP", "hR", "tP", "tI", "oP", "oS", "oI", "oF", "mP", "mS", "aP"],
            help="crystal system")

    parser.add_argument("-o", "--kpath-file", type=str, default="kpath-from-cms-2010.txt",
            help="the output kpath file")


    args = parser.parse_args()
    #
    if args.crystal_system == "cP":
        kpath = kpath_from_cms_2010(system="cP")
    elif args.crystal_system == "cI":
        kpath = kpath_from_cms_2010(system="cI")
    elif args.crystal_system == "cF":
        kpath = kpath_from_cms_2010(system="cF")

    with open(args.kpath_file, 'w') as fout:
        fout.write("%d\n" % len(kpath))
        for kpoint in kpath:
            if kpoint[4] !=  "|":
                fout.write("%f %f %f #%s %d\n" % (kpoint[0], kpoint[1], kpoint[2], kpoint[3], kpoint[4]))
            else:
                fout.write("%f %f %f #%s |\n" % (kpoint[0], kpoint[1], kpoint[2], kpoint[3]))

if __name__ == "__main__":
    main()
