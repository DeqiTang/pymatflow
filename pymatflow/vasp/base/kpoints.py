#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np

class vasp_kpoints:
    """
    """
    def __init__(self):
        self.kpoints_mp = [1, 1, 1, 0, 0, 0]

    def to_kpoints(self, fout):
        if self.kpoints_option == "automatic":
            fout.write("K-POINTS automatic\n")
            fout.write("0\n") # automatically generate the grid
            fout.write("Gamma\n") # Gamma centered
            fout.write("%d %d %d\n" % (self.kpoints_mp[0], self.kpoints_mp[1], self.kpoints_mp[2]))
            fout.write("%d %d %d\n" % (self.kpoints_mp[3], self.kpoints_mp[4], self.kpoints_mp[5])) # usually set them to 0
        elif self.kpoints_option == "bands":
            # use self.kpath
            fout.write("K-POINTS for band structure\n")
            fout.write("%d\n" % self.kpath_intersections)
            fout.write("Line-mode\n")
            fout.write("rec\n") # in reciprocal coordinates
            for i in range(len(self.kpath) - 1):
                if self.kpath[i][4] != "|" and type(self.kpath[i][4] == int):
                    fout.write(" %f %f %f !%s\n" % (self.kpath[i][0], self.kpath[i][1], self.kpath[i][2], self.kpath[i][3]))
                    fout.write(" %f %f %f !%s\n" % (self.kpath[i+1][0], self.kpath[i+1][1], self.kpath[i+1][2], self.kpath[i+1][3]))
                    fout.write("\n")
        else:
            pass


    def to_string(self):
        out = ""
        if self.kpoints_option == "automatic":
            out += "K-POINTS automatic\n"
            out += "0\n" # automatically generate the grid
            out += "Gamma\n" # Gamma centered
            out += "%d %d %d\n" % (self.kpoints_mp[0], self.kpoints_mp[1], self.kpoints_mp[2])
            out += "%d %d %d\n" % (self.kpoints_mp[3], self.kpoints_mp[4], self.kpoints_mp[5]) # usually set them to 0
        elif self.kpoints_option == "bands":
            # use self.kpath
            out += "K-POINTS for band structure\n"
            out += "%d\n" % self.kpath_intersections
            out += "Line-mode\n"
            out += "rec\n" # in reciprocal coordinates
            for i in range(len(self.kpath) - 1):
                if self.kpath[i][4] != "|" and type(self.kpath[i][4] == int):
                    out += " %f %f %f !%s\n" % (self.kpath[i][0], self.kpath[i][1], self.kpath[i][2], self.kpath[i][3])
                    out += " %f %f %f !%s\n" % (self.kpath[i+1][0], self.kpath[i+1][1], self.kpath[i+1][2], self.kpath[i+1][3])
                    out += "\n"
        else:
            pass
        return out

    def set_kpoints(self, kpoints_mp=[1, 1, 1, 0, 0, 0], option="automatic", kpath=None, kpath_intersections=15):
        """
        option: automatic, bands

        kpath:
            the high symmetry k point path used in bands structure calculation
            in format like this:

            [[kx, ky, kz, label, connect_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', 15], ...]

            if connect_indicator in a kpoint is an integer, then it will connect to the following point
            through the number of kpoints defined by connect_indicator.

            if connect_indicator in a kpoint is '|', then it will not connect to the following point,
        """
        self.kpoints_option = option
        self.kpoints_mp = kpoints_mp
        self.kpath = kpath
        self.kpath_intersections = kpath_intersections
