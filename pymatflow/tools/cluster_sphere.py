#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
from emuhelper.base.cluster import cluster_xyz

"""
usage:
    cluster_sphere input.xyz output.xyz radius
    radius in unit of angstrom
"""

def main():
    xyz = cluster_xyz(sys.argv[1])
    radius = float(sys.argv[3])
    xyz.build_cluster_sphere(radius)
    xyz.to_xyz(sys.argv[2])




if __name__ == '__main__':
    main()
