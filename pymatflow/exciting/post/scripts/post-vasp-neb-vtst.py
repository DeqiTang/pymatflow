#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of relax running", type=str, default="matflow-running")

    args = parser.parse_args()


    # making traj: initial and current
    os.chdir(args.directory)
    dirs_int = []
    for i in os.listdir():
        if i.isdecimal() == True:
            dirs_int.append(int(i))
    total_images = max(dirs_int) + 1 # including initial and final
    
    for i in range(total_images):
        if i == 0 or i == (total_images-1):
            os.system("sflow convert -i %.2d/POSCAR -o %.2d/POSCAR.xyz" % (i, i))
        else:
            os.system("sflow convert -i %.2d/POSCAR -o %.2d/POSCAR.xyz" % (i, i))
            os.system("sflow convert -i %.2d/CONTCAR -o %.2d/CONTCAR.xyz" % (i, i))
    for i in range(total_images):
        # initial traj
        if i == 0:
            os.system("cat %.2d/POSCAR.xyz > traj-initial.xyz" % (i))
        else:
            os.system("cat %.2d/POSCAR.xyz >> traj-initial.xyz" % (i))        
        # current traj
        if i == 0:
            os.system("cat %.2d/POSCAR.xyz > traj-current.xyz" % (i))        
        elif i == (total_images-1):
            os.system("cat %.2d/POSCAR.xyz >> traj-current.xyz" % (i))
        else:
            os.system("cat %.2d/CONTCAR.xyz >> traj-current.xyz" % (i))
    # end making the traj
    
    