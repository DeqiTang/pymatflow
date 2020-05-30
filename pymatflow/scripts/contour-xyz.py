#!/usr/bin/env python

import os
import sys
import argparse

from pymatflow.cmd.structflow import read_structure

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    parser.add_argument("-o", "--output", type=str, default="./contour",
        help="prefix of the output image file name")
    
    parser.add_argument("--atoms", nargs='+', type=int,
        help="list of fixed atoms, index start from 1, default is all atoms")

    parser.add_argument("--around-z", type=float, nargs=3,
        help="select atoms around specified z in Angstrom with tolerance, like this --around-z 10 -0.5 0.5")

    parser.add_argument("--dgrid3d", type=int, nargs=3,
        default=[100, 100, 4],
        help="used by gnuplot to set dgrid3d int, int, int")

    parser.add_argument("--levels", type=int, default=10,
        help="levels of the color map or color bar")

    parser.add_argument("--ngridx", type=int, default=200,
        help="ngridx to plot the contour for irregularly spaced data")

    parser.add_argument("--ngridy", type=int, default=200,
        help="ngridy to plot the contour for irregularly spaced data")        

    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()


    structure = read_structure(args.input)

    if args.atoms == None and args.around_z == None:
        atoms_index_from_1 = range(1, len(structure.atoms)+1)
    elif args.atoms == None and args.around_z != None:
        atoms_index_from_1 = []
        for i in range(len(structure.atoms)):
            if structure.atoms[i].z > (args.around_z[0] + args.around_z[1]) and structure.atoms[i].z < (args.around_z[0] + args.around_z[2]):
                atoms_index_from_1.append(i+1)
    elif args.atoms != None and args.around_z == None:
        atoms_index_from_1 = args.atoms
    elif args.atoms != None and args.around_z != None:
        print("======================================================\n")
        print("                      Warning\n")
        print("--atoms and --around-z can not be used at the same time.")
        sys.exit(1)
    
    data = []
    for i in atoms_index_from_1:
        #X.append(structure.atoms[i-1].x)
        #Y.append(structure.atoms[i-1].y)
        #Z.append(structure.atoms[i-1].z)
        data.append([structure.atoms[i-1].x, structure.atoms[i-1].y, structure.atoms[i-1].z])

    with open(args.output+".data", "w") as fout:
        fout.write("# format: x y z\n")
        for d in data:
            fout.write("%f\t%f\t%f\n" % (d[0], d[1], d[2]))
    

    with open("plot.gnuplot", "w") as fout:
        fout.write("set term gif\n")
        fout.write("set dgrid3d %d, %d, %d\n" % (args.dgrid3d[0], args.dgrid3d[1], args.dgrid3d[2]))
        fout.write("set contour base\n")
        fout.write("unset key\n")
        fout.write("set autoscale\n")
        fout.write("set xlabel \"X\"\n")
        fout.write("set ylabel \"Y\"\n")
        fout.write("set zlabel \"Z\"\n")
        #fout.write("set output \"%s\"\n" % (args.output+".gif"))
        #fout.write("set multiplot layout 1, 2\n")

        #fout.write("set origin 0, 0\n")
        fout.write("\n\n")
        fout.write("set surface\n")
        fout.write("set pm3d hidden3d\n")
        fout.write("set output \"%s\"\n" % (args.output+".surf"+".gif"))
        fout.write("splot '%s' u 1:2:3 w pm3d\n" % (args.output+".data"))

        #fout.write("set origin 0.5, 0\n")
        fout.write("\n\n")
        fout.write("set pm3d map\n")
        fout.write("set surface\n")
        fout.write("set output \"%s\"\n" % (args.output+".map"+".gif"))
        fout.write("splot '%s' u 1:2:3 w pm3d\n" % (args.output+".data"))
    os.system("gnuplot plot.gnuplot")

    # ----------------
    # 2D contour plot
    #-----------------
    # Contour plot of irregularly spaced data
    # data is irregularly spaced data
    # we can refer to https://matplotlib.org/3.2.1/gallery/images_contours_and_fields/irregulardatagrid.html
    # to see how to build contour plot of such kind of data
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri
    
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
    x = data_np[:, 0]
    y = data_np[:, 1]
    z = data_np[:, 2]

    fig, (ax1, ax2) = plt.subplots(nrows=2)

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

    ax1.contour(xi, yi, zi, levels=args.levels, linewidths=0.5, colors='k')
    cntr1 = ax1.contourf(xi, yi, zi, levels=args.levels, cmap="RdBu_r")

    fig.colorbar(cntr1, ax=ax1)
    ax1.plot(x, y, 'ko', ms=3)
    #ax1.set(xlim=(-2, 2), ylim=(-2, 2))
    #ax1.set_title('grid and contour (%d points, %d grid points)' %
                #(npts, ngridx * ngridy))

    ax1.axis('equal') # set x y axis equal spaced

    # ----------
    # Tricontour
    # ----------
    # Directly supply the unordered, irregularly spaced coordinates
    # to tricontour.

    ax2.tricontour(x, y, z, levels=args.levels, linewidths=0.5, colors='k')
    cntr2 = ax2.tricontourf(x, y, z, levels=args.levels, cmap="RdBu_r")

    fig.colorbar(cntr2, ax=ax2)
    ax2.plot(x, y, 'ko', ms=3)
    #ax2.set(xlim=(-2, 2), ylim=(-2, 2))
    #ax2.set_title('tricontour (%d points)' % npts)
    
    ax2.axis('equal') # set x y axis equal spaced

    plt.subplots_adjust(hspace=0.5)
    plt.autoscale()
    plt.show() 
    fig.savefig(args.output+".2d-contour.png")

if __name__ == "__main__":
    main()