#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse


def main():
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
    
    os.system("mkdir -p post-processing")

    for i in range(total_images):
        # initial traj
        if i == 0:
            os.system("cat %.2d/POSCAR.xyz > post-processing/traj-initial.xyz" % (i))
        else:
            os.system("cat %.2d/POSCAR.xyz >> post-processing/traj-initial.xyz" % (i))        
        # current traj
        if i == 0:
            os.system("cat %.2d/POSCAR.xyz > post-processing/traj-current.xyz" % (i))        
        elif i == (total_images-1):
            os.system("cat %.2d/POSCAR.xyz >> post-processing/traj-current.xyz" % (i))
        else:
            os.system("cat %.2d/CONTCAR.xyz >> post-processing/traj-current.xyz" % (i))
    # end making the traj


    # build the energy barrier plot
    import numpy as np 
    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d
    
    os.system("nebresults.pl")
    
    neb_dat = np.loadtxt("neb.dat")
    
    # style: linear
    plt.plot(neb_dat[:, 1], neb_dat[:, 2], marker="o")
    plt.xlabel("Reaction Coordinate (Angstrom)")
    plt.ylabel("Energy (eV)")
    plt.savefig("post-processing/mep-style-linear.png")
    plt.close()

    
    # style: cubic
    fun_dat = interp1d(neb_dat[:, 1], neb_dat[:, 2], kind="cubic")
    denser_x = np.linspace(neb_dat[:, 1].min(), neb_dat[:, 1].max(), 30*len(neb_dat))
    plt.plot(denser_x, fun_dat(denser_x))
    plt.scatter(neb_dat[:, 1], neb_dat[:, 2], marker="o")
    plt.xlabel("Reaction Coordinate (Angstrom)")
    plt.ylabel("Energy (eV)")
    plt.savefig("post-processing/mep-style-cubic-order.png")
    plt.close()    
        
    # style: order of 5
    fun_dat = interp1d(neb_dat[:, 1], neb_dat[:, 2], kind=5)
    denser_x = np.linspace(neb_dat[:, 1].min(), neb_dat[:, 1].max(), 30*len(neb_dat))
    plt.plot(denser_x, fun_dat(denser_x))
    plt.scatter(neb_dat[:, 1], neb_dat[:, 2], marker="o")
    plt.xlabel("Reaction Coordinate (Angstrom)")
    plt.ylabel("Energy (eV)")
    plt.savefig("post-processing/mep-style-5th-order.png")
    plt.close()    
    
    

if __name__ == "__main__":
    main()