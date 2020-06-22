#!/usr/bin/env python

import os
import sys
import copy
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

from pymatflow.cmd.structflow import read_structure
from pymatflow.cmd.structflow import write_structure

from pymatflow.vasp.post.pdos import post_pdos

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str, required=True,
        help="input PARCHG file")

    parser.add_argument("--output-structure", type=str, default="parchg.cif",
        help="output stucture contained in PARCHG")

    parser.add_argument("-o", "--output", type=str, default="./stm",
        help="prefix of the output image file name")

    parser.add_argument("--dgrid3d", type=int, nargs=3,
        default=[100, 100, 4],
        help="used by gnuplot to set dgrid3d int, int, int")
    
    parser.add_argument("--levels", type=int, default=10,
        help="levels of the color map or color bar")


    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()
    
    parchg_filepath = args.input
    
    with open(parchg_filepath, "r") as fin:
        parchg = fin.readlines()
        
    for i in range(len(parchg)):
        if len(parchg[i].split()) == 0:
            first_blank_line = i
            break
    os.system("mkdir -p /tmp/pymatflow/")
    with open("/tmp/pymatflow/POSCAR", "w") as fout:
        for i in range(first_blank_line):
            fout.write(parchg[i])
    structure = read_structure("/tmp/pymatflow/POSCAR")
    write_structure(structure=structure, filepath=args.output_structure)
    
    ngxf = int(parchg[first_blank_line+1].split()[0])
    ngyf = int(parchg[first_blank_line+1].split()[1])
    ngzf = int(parchg[first_blank_line+1].split()[2])    
    
    data = np.loadtxt(parchg[first_blank_line+2:])
    data = data.reshape(ngzf, ngyf, ngxf)
    
    
    # ----------------
    # gray scale image
    # ----------------
    
    os.system("mkdir -p tmp-stm-images")
    for i in range(data.shape[0]):
        img = data[i, ::-1, ::]
        img = (img-img.min()) / (img.max() - img.min()) * 255
        # need to do a transform when the cell is not Orthorhombic
        # skew the image
        ax = plt.axes() #plt.figure()
        im = ax.imshow(img, cmap="gray")
        #im = plt.imshow(data[i, :, :], cmap="gray")
        a = np.array(structure.cell[0])
        b = np.array(structure.cell[1])
        cosangle = a.dot(b)/(np.linalg.norm(a) * np.linalg.norm(b))
        angle = np.arccos(cosangle) * 180 / np.pi
        trans_data = ax.transData + mtransforms.Affine2D().skew_deg(90-angle, 0)
        im.set_transform(trans_data)
        plt.colorbar(im)
        ax.autoscale()
        plt.tight_layout()
        plt.savefig("tmp-stm-images/"+args.output+"%d.gray.png" % i)
        plt.close()
        
    # ----------------
    # 2D contour plot
    #-----------------
    
    nx = np.linspace(0, 1, ngxf)
    ny = np.linspace(0, 1, ngyf)
    X, Y = np.meshgrid(nx, ny) # now this Mesh grid cannot be used directly, we have to calc the real x y for it
    for xi in range(len(nx)):
        for yi in range(len(ny)):
            X[yi, xi] = structure.cell[0][0] * nx[xi] + structure.cell[1][0] * ny[yi]
            Y[yi, xi] = structure.cell[0][1] * nx[xi] + structure.cell[1][1] * ny[yi]
    
    Z = data[0, :, :]
    Z = (Z-Z.min()) / (Z.max() - Z.min()) * 255
    # fill color, three color are divided into three layer(6)
    # cmap = plt.cm.hot means using thermostat plot(graduated red yellow)
    #cset = plt.contourf(X, Y, Z, levels=args.levels, cmap=plt.cm.hot)
    cset = plt.contourf(X, Y, Z, levels=args.levels, cmap=plt.cm.gray)
    contour = plt.contour(X, Y, Z, levels=[20, 40], colors='k')
    plt.colorbar(cset)
    plt.autoscale()
    plt.tight_layout()
    #plt.show()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(args.output+".2d-contour.png")
    plt.close()
        
    

if __name__ == "__main__":
    main()