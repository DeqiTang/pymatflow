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
    
    parser.add_argument("-i", "--input", type=str, nargs=3, required=True,
        help="input vasp *CHG* file, -i TOTAL PART1 PART2")

    parser.add_argument("--output-structure", type=str, default="diff-chg",
        help="output stucture contained in *CHG*")

    parser.add_argument("-o", "--output", type=str, default="diff-chg",
        help="prefix of the output image file name")
    
    parser.add_argument("--levels", type=int, default=10,
        help="levels of the color map or color bar")

    parser.add_argument("-z", "--z", type=float, default=1,
        help="a value between 0 and 1, indicat height of in z direction to print the plot")
        
    parser.add_argument("--cmap", type=str, default="gray",
        choices=["gray", "hot", "afmhot", "Spectral", "plasma", "magma", "hsv", "rainbow", "brg"])
        
    parser.add_argument("--abscissa", type=str, nargs="+", default=["a", "b", "c"], 
        choices=["a", "b", "c"], 
        help="choose the direction to do the dimension reduction")
        
    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()
    
    chg_filepath = args.input
    
    chg = []
    for i in range(3):
        with open(chg_filepath[i], "r") as fin:
            chg.append(fin.readlines())
    
    first_blank_line = []
    for i in range(3):
        for j in range(len(chg[i])):
            if len(chg[i][j].split()) == 0:
                first_blank_line.append(j)
                break
    
    first_augmentation_line = [None, None, None]
    for i in range(3):
        for j in range(len(chg[i])):
            if "augmentation" in chg[i][j]:
                first_augmentation_line[i] = j
                break
                
    os.system("mkdir -p /tmp/pymatflow/")
    for i in range(3):
        with open("/tmp/pymatflow/POSCAR", "w") as fout:
            for j in range(first_blank_line[i]):
                fout.write(chg[i][j])
            
        structure = read_structure("/tmp/pymatflow/POSCAR")
        if i == 0:
            write_structure(structure=structure, filepath=args.output_structure+".total.cif")
        elif i == 1:
            write_structure(structure=structure, filepath=args.output_structure+".part1.cif")
        elif i == 2:
            write_structure(structure=structure, filepath=args.output_structure+".part2.cif")
    
    a = np.linalg.norm(structure.cell[0])
    b = np.linalg.norm(structure.cell[1])
    c = np.linalg.norm(structure.cell[2])
    
    # assume three *CHG* have the same ngxf and ngyf ngzf
    ngxf = int(chg[0][first_blank_line[0]+1].split()[0])
    ngyf = int(chg[0][first_blank_line[0]+1].split()[1])
    ngzf = int(chg[0][first_blank_line[0]+1].split()[2])    
    
    data_iii = []
    for i in range(3):
        if first_augmentation_line[i] == None:
            tmp_str = "".join(chg[i][first_blank_line[i]+2:])
            data_iii.append(np.fromstring(tmp_str, sep="\n").reshape(ngzf, ngyf, ngxf))
        else:
            tmp_str = "".join(chg[i][first_blank_line[i]+2:first_augmentation_line[i]])
            data_iii.append(np.fromstring(tmp_str, sep="\n").reshape(ngzf, ngyf, ngxf))
        #data = data.reshape(ngzf, ngyf, ngxf)
    
    
    # data_iii dimension reduction
    data_sub = np.array(data_iii[0]) - np.array(data_iii[1]) - np.array(data_iii[2])
    # the unit of value is actually not physical now!
    cell_volume = np.dot(np.cross(np.array(structure.cell[0]), np.array(structure.cell[1])), np.array(structure.cell[2]))
    cell_volume_per_unit = cell_volume / (ngzf * ngyf * ngxf)
    
    # value in Vasp *CHG* are \rho(r)_of_electrons * Volume_of_cell, so we should divide it by cell_volume here and time it with cell_volume_per_unit
    # to get the number of electrons per divided unit
    total_electrons = np.sum(data_iii[0]) / cell_volume * cell_volume_per_unit
    #
    
    print("======================================================\n")
    print("           Information collected\n")
    print("------------------------------------------------------\n")
    print("cell volume: %f (A^3)\n" % cell_volume)
    print("total electrons: %f\n" % total_electrons)
    
    
    
    # unit of data_red_? is e/Anstrom, namely number of electrons per Angstrom
    data_red_a = []
    data_red_b = []
    data_red_c = []
    if "c" in args.abscissa:
        factor = cell_volume_per_unit / cell_volume
        len_ci = c / ngzf
        for ci in range(data_sub.shape[0]):
            tmp = 0
            for bi in range(data_sub.shape[1]):
                tmp += np.sum(data_sub[ci, bi, :])
            nelect_ci = tmp * factor
            rho_line = nelect_ci / len_ci
            data_red_c.append(rho_line)
    if "b" in args.abscissa:
        factor = cell_volume_per_unit / cell_volume
        len_bi = b / ngyf    
        for bi in range(data_sub.shape[1]):
            tmp = 0
            for ai in range(data_sub.shape[2]):
                tmp += np.sum(data_sub[:, bi, ai])
            nelect_bi = tmp * factor
            rho_line = nelect_bi / len_bi                
            data_red_b.append(rho_line)
    if "a" in args.abscissa:
        factor = cell_volume_per_unit / cell_volume
        len_ai = a / ngxf        
        for ai in range(data_sub.shape[2]):
            tmp = 0
            for ci in range(data_sub.shape[0]):
                tmp += np.sum(data_sub[ci, :, ai])
            nelect_ai = tmp * factor
            rho_line = nelect_ai / len_ai                   
            data_red_a.append(rho_line)    

    # output the data and make the plot
    if "c" in args.abscissa:
        with open(args.output+".1d.c.data", 'w') as fout:
            fout.write("#c(angstrom) \Delta(rho) (number of electron per Angstrom)\n")
            c_coord = np.linspace(0, c, len(data_red_c))
            for i in range(len(data_red_c)):
                fout.write("%f %f\n" % (c_coord[i], data_red_c[i]))
        plt.plot(np.linspace(0, c, len(data_red_c)), data_red_c)                
        plt.ylabel(r"$\Delta\rho (e/\AA)$")
        plt.tight_layout()
        plt.savefig(args.output+".1d.c.png")
        plt.close()                
    if "b" in args.abscissa:
        with open(args.output+".1d.b.data", 'w') as fout:
            fout.write("#b(angstrom) \Delta(rho) (number of electron per Angstrom)\n")
            b_coord = np.linspace(0, b, len(data_red_b))
            for i in range(len(data_red_b)):
                fout.write("%f %f\n" % (b_coord[i], data_red_b[i]))        
        plt.plot(np.linspace(0, b, len(data_red_b)), data_red_b)    
        plt.ylabel(r"$\Delta\rho (e/\AA)$")              
        plt.tight_layout()
        plt.savefig(args.output+".1d.b.png")
        plt.close()                
    if "a" in args.abscissa:
        with open(args.output+".1d.a.data", 'w') as fout:
            fout.write("#a(angstrom) \Delta(rho) (number of electron per Angstrom)\n")
            a_coord = np.linspace(0, a, len(data_red_a))
            for i in range(len(data_red_a)):
                fout.write("%f %f\n" % (a_coord[i], data_red_a[i]))
        plt.plot(np.linspace(0, a, len(data_red_a)), data_red_a)                
        plt.ylabel(r"$\Delta\rho (e/\AA)$")            
        plt.tight_layout()
        plt.savefig(args.output+".1d.a.png")
        plt.close()
    
    
    # image data for subtracted data
    # -------------------------------------------------------
    # gray scale image only for z direction
    # may not work for triclinic and monoclinic crystal system
    # -------------------------------------------------------
    
    
    zi = int((data_sub.shape[0]-1) * args.z)
    #img = data_sub[i, ::-1, ::]
    img = data_sub[zi, ::, ::]
    img = (img-img.min()) / (img.max() - img.min()) * 255
    # need to do a transform when the cell is not Orthorhombic
    # skew the image
    a = np.array(structure.cell[0])
    b = np.array(structure.cell[1])
    cosangle = a.dot(b)/(np.linalg.norm(a) * np.linalg.norm(b))
    angle = np.arccos(cosangle) * 180 / np.pi        
    ax = plt.axes() #plt.figure()
    n1 = ngxf
    n2 = ngyf
    n1_right = n1
    n1_left = -(n2 * np.tan((angle - 90) / 180 * np.pi))
    #im = ax.imshow(img, cmap="gray", extent=[n1_left, n1_right, 0, n2], interpolation="none", origin="lower", clip_on=True)
    im = ax.imshow(img, cmap="gray", extent=[0, n1, 0, n2], interpolation="none", origin="lower", clip_on=True)
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
    plt.savefig(args.output+"-z-%f.grayscale.plot.png" % args.z)
    #plt.imsave(args.output+"-z-%f.grayscale.image.png" % args.z, arr=img, cmap="gray") # img is not Affine transoformed here!!!
    plt.close()
        
    # -----------------------------------------------------------------------------
    # 2D contour plot
    #------------------------------------------------------------------------------
    
    nx = np.linspace(0, 1, ngxf)
    ny = np.linspace(0, 1, ngyf)
    X, Y = np.meshgrid(nx, ny) # now this Mesh grid cannot be used directly, we have to calc the real x y for it
    for xi in range(len(nx)):
        for yi in range(len(ny)):
            X[yi, xi] = structure.cell[0][0] * nx[xi] + structure.cell[1][0] * ny[yi]
            Y[yi, xi] = structure.cell[0][1] * nx[xi] + structure.cell[1][1] * ny[yi]
    
    Z = data_sub[zi, :, :]
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
    plt.axis("equal") # set axis equally spaced
    #plt.show()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(args.output+".2d-contour-z-%f.png" % args.z)
    plt.close()
        
    

if __name__ == "__main__":
    main()