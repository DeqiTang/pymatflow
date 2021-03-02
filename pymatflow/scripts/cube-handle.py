#!/usr/bin/env python

import os
import sys
import copy
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

from pymatflow.base.element import element as elem
from pymatflow.structure.crystal import crystal
from pymatflow.cmd.structflow import read_structure
from pymatflow.cmd.structflow import write_structure

from pymatflow.vasp.post.pdos import post_pdos

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str, required=True,
        help="input cube file")

    parser.add_argument("--output-structure", type=str, default="cube.cif",
        help="output stucture contained in PARCHG")

    parser.add_argument("-o", "--output", type=str, default="cube",
        help="prefix of the output image file name")
    
    parser.add_argument("--levels", type=int, default=10,
        help="levels of the color map or color bar")

    parser.add_argument("-z", "--z", type=float, default=1,
        help="a value between 0 and 1, indicat height of in z direction to print the plot")
        
    parser.add_argument("--cmap", type=str, default="gray",
        choices=["gray", "hot", "afmhot", "Spectral", "plasma", "magma", "hsv", "rainbow", "brg"])
    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()
    
    cube_filepath = args.input
    
    with open(cube_filepath, "r") as fin:
        cube = fin.readlines()
        
    ngridx = int(cube[3].split()[0])
    ngridy = int(cube[4].split()[0])
    ngridz = int(cube[5].split()[0])            
        
    # read structure info
    natom = abs(int(cube[2].split()[0])) # it might be negative, if MO infor are included in cube file
    bohr_to_angstrom = 0.529177249
    structure = crystal()
    structure.cell = []
    for i in range(3):
        tmp = []
        for j in range(3):
            tmp.append(int(cube[i+3].split()[0]) * float(cube[i+3].split()[j+1]) * bohr_to_angstrom)
        structure.cell.append(tmp)
    atoms_list = []
    for i in range(natom):
        atomic_number = int(cube[i+6].split()[0])
        for e in elem:
            if elem[e].number == atomic_number:
                label = e
        atoms_list.append([
            label,
            float(cube[i+6].split()[2]) * bohr_to_angstrom,
            float(cube[i+6].split()[3]) * bohr_to_angstrom,
            float(cube[i+6].split()[4]) * bohr_to_angstrom,
        ])
    structure.get_atoms(atoms_list)
    # end read structure info
    write_structure(structure=structure, filepath=args.output_structure)


    # read grid value
    tmp_str = "".join(cube[natom+6:])
    data = np.fromstring(tmp_str, sep="\n")
        
    #data = data.reshape(ngridz, ngridy, ngridx)
    # grid data in cube is iterated in different compared to *CHG* of vasp
    data = data.reshape(ngridx, ngridy, ngridz) 
    # charge data in cube file is in shape (ngridx, ngridy, ngridz)
    # while charge in *CHG* file is in shape (ngzf, ngyf, ngxf)
    # they are different!
    
    
    # -------------------------------------------------------
    # matrix image only for z direction
    # may not work for triclinic and monoclinic crystal system
    # -------------------------------------------------------
    
    
    zi = int((data.shape[2]-1) * args.z)
    img = data[::, ::, zi].T #.reshape(ngridy, ngridx) should be transpose here but not reshape
    img = (img-img.min()) / (img.max() - img.min()) * 255
    # need to do a transform when the cell is not Orthorhombic
    # skew the image
    a = np.array(structure.cell[0])
    b = np.array(structure.cell[1])
    cosangle = a.dot(b)/(np.linalg.norm(a) * np.linalg.norm(b))
    angle = np.arccos(cosangle) * 180 / np.pi        
    ax = plt.axes() #plt.figure()
    n1 = ngridx
    n2 = ngridy
    n1_right = n1
    n1_left = -(n2 * np.tan((angle - 90) / 180 * np.pi))
    #im = ax.imshow(img, cmap="gray", extent=[0, n1, 0, n2], interpolation="none", origin="lower", clip_on=True)
    im = ax.imshow(img, cmap=args.cmap, extent=[0, n1, 0, n2], interpolation="none", origin="lower", clip_on=True)
    #im = plt.imshow(data[i, :, :], cmap="gray")
    trans_data = mtransforms.Affine2D().skew_deg(90-angle, 0) + ax.transData
    im.set_transform(trans_data)
    # display intended extent of the image
    x1, x2, y1, y2 = im.get_extent()
    # do not view the line, but it is needed to be plot so the intended image is dispalyed completely
    ax.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], "y--", transform=trans_data, visible=False) 
    #ax.set_xlim(n1_left, n1_right)
    #ax.set_ylim(0, n2)
    plt.colorbar(im)
    ax.autoscale()
    plt.tight_layout()
    plt.savefig(args.output+"-z-%f.matrix-image.png" % args.z)
    plt.close()
        
    # -----------------------------------------------------------------------------
    # 2D contour plot
    #------------------------------------------------------------------------------
    
    nx = np.linspace(0, 1, ngridx)
    ny = np.linspace(0, 1, ngridy)
    X, Y = np.meshgrid(nx, ny) # now this Mesh grid cannot be used directly, we have to calc the real x y for it
    for xi in range(len(nx)):
        for yi in range(len(ny)):
            X[yi, xi] = structure.cell[0][0] * nx[xi] + structure.cell[1][0] * ny[yi]
            Y[yi, xi] = structure.cell[0][1] * nx[xi] + structure.cell[1][1] * ny[yi]
    
    Z = data[::, ::, zi].T #.reshape(ngridy, ngridx) should be transpose here but not reshape
    Z = (Z-Z.min()) / (Z.max() - Z.min()) * 255
    # fill color, three color are divided into three layer(6)
    # cmap = plt.cm.hot means using thermostat plot(graduated red yellow)
    #cset = plt.contourf(X, Y, Z, levels=args.levels, cmap=plt.cm.hot)
    #cset = plt.contourf(X, Y, Z, levels=args.levels, cmap=plt.cm.gray)
    cset = plt.contourf(X, Y, Z, levels=args.levels, cmap=args.cmap)
    contour = plt.contour(X, Y, Z, levels=[20, 40], colors='k')
    plt.colorbar(cset)
    plt.autoscale()
    plt.tight_layout()
    #plt.show()
    plt.axis("equal") # set axis equally spaced
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(args.output+".2d-contour-z-%f.png" % args.z)
    plt.close()
        
    

if __name__ == "__main__":
    main()