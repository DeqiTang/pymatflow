#!/usr/bin/env python

import os
import sys
import copy
import argparse

from pymatflow.cmd.structflow import read_structure

from pymatflow.vasp.post.pdos import post_pdos

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str, nargs="+", required=True,
            help="input structure file and vasprun.xml, you can specify two vasprun.xml(one for scf and one for nscf). eg. -i structure.cif vasprun.xml.scf vasprun.xml.nscf")

    parser.add_argument("-o", "--output", type=str, default="./contour-ldos",
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

    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()


    structure = read_structure(args.input[0])

    pdos = post_pdos()
    if len(args.input) == 2:
        efermi_from = "nscf"
        pdos.get_vasprun(args.input[1])
        pdos.get_efermi(vasprun=args.input[1])
    elif len(args.input) == 3:
        efermi_from = "scf"
        pdos.get_vasprun(args.input[2])
        pdos.get_efermi(vasprun=args.input[1])
    #pdos.export(directory=args.directory, engine=args.engine, plotrange=args.plotrange)



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
    
    # data: [x, y, z, ldos]
    data = []
    for i in atoms_index_from_1:
        ldos = 0.0
        #if pdos.magnetic_status == "non-soc-ispin-1":
        #    pdos.data[i-1]
        #elif pdos.magnetic_status == "non-soc-ispin-2":
        #    pdos.data[i-1]
        #elif pdos.magnetic_status == "soc-ispin-1" or pdos.magnetic_status == "soc-ispin-2":
        #    pdos.data[i-1]
        for item in pdos.data[i-1]:
            if item == "ion":
                continue
            for itm in pdos.data[i-1][item]:
                if itm == "energy":
                    continue
                ldos += sum(pdos.data[i-1][item][itm])
        data.append([structure.atoms[i-1].x, structure.atoms[i-1].y, structure.atoms[i-1].z, ldos])


    with open(args.output+".data", "w") as fout:
        fout.write("# format: x y z ldos\n")
        for d in data:
            fout.write("%f\t%f\t%f\t%f\n" % (d[0], d[1], d[2], d[3]))
    


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
    x = data_np[:, 0] # x
    y = data_np[:, 1] # y
    z = data_np[:, 3] # ldos

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

    plt.close()
    ax = plt.axes(projection='3d')    
    cset = ax.plot_trisurf(x, y, z, cmap='rainbow')
    plt.colorbar(cset)
    plt.autoscale()
    plt.tight_layout()
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('ldos')    
    plt.savefig(args.output+".3d-trisurf.png")

if __name__ == "__main__":
    main()