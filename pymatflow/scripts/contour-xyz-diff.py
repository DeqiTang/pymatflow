#!/usr/bin/env python

import os
import sys
import copy
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
    
from pymatflow.cmd.structflow import read_structure

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str, nargs=2, required=True,
            help="input initial and final structure file, eg. -i initial.cif final.cif")

    parser.add_argument("-o", "--output", type=str, default="./contour-diff",
        help="prefix of the output image file name")
    
    parser.add_argument("--atoms", nargs='+', type=int,
        help="list of fixed atoms, index start from 1, default is all atoms")

    parser.add_argument("--around-z", type=float, nargs=3,
        help="select atoms around specified z in Angstrom with tolerance, like this --around-z 10 -0.5 0.5")


    parser.add_argument("--levels", type=int, default=10,
        help="levels of the color map or color bar")

    parser.add_argument("--ngridx", type=int, default=200,
        help="ngridx to plot the contour for irregularly spaced data")

    parser.add_argument("--ngridy", type=int, default=200,
        help="ngridy to plot the contour for irregularly spaced data")        

    parser.add_argument("--diff", type=str, default="z",
        choices=["x", "y", "z", "xyz"],
        help="choose to plot the diff of x or y or z, or xyz which means displacement")
    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()


    initial_structure = read_structure(args.input[0])
    final_structure = read_structure(args.input[1])

    if args.atoms == None and args.around_z == None:
        atoms_index_from_1 = range(1, len(initial_structure.atoms)+1)
    elif args.atoms == None and args.around_z != None:
        atoms_index_from_1 = []
        for i in range(len(initial_structure.atoms)):
            if initial_structure.atoms[i].z > (args.around_z[0] + args.around_z[1]) and initial_structure.atoms[i].z < (args.around_z[0] + args.around_z[2]):
                atoms_index_from_1.append(i+1)
    elif args.atoms != None and args.around_z == None:
        atoms_index_from_1 = args.atoms
    elif args.atoms != None and args.around_z != None:
        print("======================================================\n")
        print("                      Warning\n")
        print("--atoms and --around-z can not be used at the same time.")
        sys.exit(1)
    
    initial_data = []
    for i in atoms_index_from_1:
        #X.append(structure.atoms[i-1].x)
        #Y.append(structure.atoms[i-1].y)
        #Z.append(structure.atoms[i-1].z)
        initial_data.append([initial_structure.atoms[i-1].x, initial_structure.atoms[i-1].y, initial_structure.atoms[i-1].z])

    final_data = []
    for i in atoms_index_from_1:
        #X.append(structure.atoms[i-1].x)
        #Y.append(structure.atoms[i-1].y)
        #Z.append(structure.atoms[i-1].z)
        final_data.append([final_structure.atoms[i-1].x, final_structure.atoms[i-1].y, final_structure.atoms[i-1].z])


    # data: [x, y, z, diff(x|y|z|xyz)]
    data = copy.deepcopy(final_data)
    for i in range(len(final_data)):
        if args.diff == "x":
            diff = final_data[i][0] - initial_data[i][0]
        elif args.diff == "y":
            diff = final_data[i][1] - initial_data[i][1]   
        elif args.diff == "z":
            diff = final_data[i][2] - initial_data[i][2]
        elif args.diff == "xyz":
            diff =  np.sqrt((final_data[i][0] - initial_data[i][0])**2 + (final_data[i][1] - initial_data[i][1])**2 + (final_data[i][2] - initial_data[i][2])**2)
        else:
            pass
        data[i].append(diff)

    with open(args.output+".data", "w") as fout:
        fout.write("# format: x y z diff(%s)\n" % (args.diff))
        for d in data:
            fout.write("%f\t%f\t%f\t%f\n" % (d[0], d[1], d[2], d[3]))
    


    # ----------------
    # 2D contour plot
    #-----------------
    # Contour plot of irregularly spaced data
    # data is irregularly spaced data
    # we can refer to https://matplotlib.org/3.2.1/gallery/images_contours_and_fields/irregulardatagrid.html
    # to see how to build contour plot of such kind of data

    
    data_np = np.array(data)

    # fill color, three color are divided into three layer(6)
    # cmap = plt.cm.hot means using thermostat plot(graduated red yellow)
    #cset = plt.contourf(X, Y, Z, levels=10, cmap=plt.cm.hot)
    #contour = plt.contour(X, Y, Z, colors='k')
    #plt.colorbar(cset)
    #plt.autoscale()
    #plt.tight_layout()
    #plt.show()
    #plt.xlabel('x')
    #plt.ylabel('y')
    #plt.savefig(args.output+".2d-contour.png")
    #plt.close()    
    
    ngridx = args.ngridx
    ngridy = args.ngridy
    x = data_np[:, 0] # x
    y = data_np[:, 1] # y
    z = data_np[:, 3] # diff(x|y|z)

    #fig, (ax1, ax2) = plt.subplots(nrows=2)
    fig = plt.figure()
    # -----------------------
    # Interpolation on a grid
    # -----------------------
    # A contour plot of irregularly spaced data coordinates
    # via interpolation on a grid.

    # Create grid values first.
    xi = np.linspace(np.min(x), np.max(x), ngridx)
    yi = np.linspace(np.min(y), np.max(y), ngridy)

    # Linearly interpolate the data (x, y) on a grid defined by (xi, yi).
    triang = tri.Triangulation(x, y)
    interpolator = tri.LinearTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    zi = interpolator(Xi, Yi)

    # Note that scipy.interpolate provides means to interpolate data on a grid
    # as well. The following would be an alternative to the four lines above:
    #from scipy.interpolate import griddata
    #zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')

    #ax1.contour(xi, yi, zi, levels=args.levels, linewidths=0.5, colors='k')

    #cntr1 = ax1.contourf(xi, yi, zi, levels=args.levels, cmap="RdBu_r")

    #fig.colorbar(cntr1, ax=ax1)
    #ax1.plot(x, y, 'ko', ms=3)
    #ax1.set(xlim=(-2, 2), ylim=(-2, 2))
    #ax1.set_title('grid and contour (%d points, %d grid points)' %
                #(npts, ngridx * ngridy))

    #ax1.axis('equal') # set x y axis equal spaced

    # ----------
    # Tricontour
    # ----------
    # Directly supply the unordered, irregularly spaced coordinates
    # to tricontour.

    #ax2.tricontour(x, y, z, levels=args.levels, linewidths=0.5, colors='k')
    plt.tricontour(x, y, z, levels=args.levels, linewidths=0.5, colors='k')
    #cntr2 = ax2.tricontourf(x, y, z, levels=args.levels, cmap="RdBu_r")
    cntr2 = plt.tricontourf(x, y, z, levels=args.levels, cmap="RdBu_r")

    #fig.colorbar(cntr2, ax=ax2)
    fig.colorbar(cntr2)
    #ax2.plot(x, y, 'ko', ms=3)
    #plt.plot(x, y, "ko", ms=3)
    #ax2.set(xlim=(-2, 2), ylim=(-2, 2))
    #ax2.set_title('tricontour (%d points)' % npts)
    
    #ax2.axis('equal') # set x y axis equal spaced
    plt.axis("equal")

    #plt.subplots_adjust(hspace=0.5)
    plt.autoscale()
    plt.show() 
    fig.savefig(args.output+".2d-contour.png")

if __name__ == "__main__":
    main()