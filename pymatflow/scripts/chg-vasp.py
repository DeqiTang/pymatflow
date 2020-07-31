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
        help="input vasp *CHG* file")

    parser.add_argument("--output-structure", type=str, default="chg.cif",
        help="output stucture contained in PARCHG")
        
    parser.add_argument("-o", "--output", type=str, default="chg",
        help="prefix of the output image file name")        
    
    parser.add_argument("--format", type=str, default="png",
        choices=["png"],
        help="format of the image")        
    
    parser.add_argument("--output-option", type=int, nargs="+",
        default=[1, 2, 3],
        help="choose to output images of many kinds! (1)->grayscale image in z direction with scale bar; (2)->2D contour plot; (3)->grayscale image without scale bar")
    
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
    
    chg_filepath = args.input
    
    with open(chg_filepath, "r") as fin:
        chg = fin.readlines()
        
    for i in range(len(chg)):
        if len(chg[i].split()) == 0:
            first_blank_line = i
            break
    
    first_augmentation_line = None
    for i in range(len(chg)):
        if "augmentation" in chg[i]:
            first_augmentation_line = i
            break
                
    os.system("mkdir -p /tmp/pymatflow/")
    with open("/tmp/pymatflow/POSCAR", "w") as fout:
        for i in range(first_blank_line):
            fout.write(chg[i])
            
    structure = read_structure("/tmp/pymatflow/POSCAR")
    write_structure(structure=structure, filepath=args.output_structure)
    
    is_orthogonal = False
    cos_ab = np.dot(np.array(structure.cell[0]), np.array(structure.cell[1]))
    cos_ac = np.dot(np.array(structure.cell[0]), np.array(structure.cell[2]))
    cos_bc = np.dot(np.array(structure.cell[1]), np.array(structure.cell[2]))
    if cos_ab == cos_ac == cos_bc == 0:
        is_orthogonal = True
    
    ngxf = int(chg[first_blank_line+1].split()[0])
    ngyf = int(chg[first_blank_line+1].split()[1])
    ngzf = int(chg[first_blank_line+1].split()[2])    
    
    if first_augmentation_line == None:
        #data = np.loadtxt(chg[first_blank_line+2:])
        tmp_str = "".join(chg[first_blank_line+2:])
        data = np.fromstring(tmp_str, sep="\n")
    else:
        #data = np.loadtxt(chg[first_blank_line+2:first_augmentation_line])
        tmp_str = "".join(chg[first_blank_line+2:first_augmentation_line])
        data = np.fromstring(tmp_str, sep="\n")
        
    data = data.reshape(ngzf, ngyf, ngxf)
    
    
    if 1 in args.output_option:
        # -------------------------------------------------------
        # gray scale image only for z direction
        # may not work for triclinic and monoclinic crystal system
        # -------------------------------------------------------
    
    
        zi = int((data.shape[0]-1) * args.z)
        #img = data[i, ::-1, ::]
        img = data[zi, ::, ::]
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
        plt.savefig(args.output+"-z-%f.grayscale.plot.%s" % (args.z, args.format))
        #plt.imsave(args.output+"-z-%f.grayscale.image.png" % args.z, arr=img, cmap="gray") # img is not Affine transoformed here!!!
        plt.close()
        
    if 2 in args.output_option:
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
    
        Z = data[zi, :, :]
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
        plt.savefig(args.output+".2d-contour-z-%f.%s" % (args.z, args.format))
        plt.close()
        
    
    if (3 in args.output_option) and is_orthogonal == True:
        # -------------------------------------------------------
        # gray scale image only for orthogonal crystal system
        # -------------------------------------------------------
        zi = int((data.shape[0]-1) * args.z)
        #img = data[i, ::-1, ::]
        img = data[zi, ::, ::]
        img = (img-img.min()) / (img.max() - img.min()) * 255
        
        import PIL.Image as Image
        
        img = Image.fromarray(img.astype('uint8'))
        img = img.convert('L') # gryscale
        img.save(args.output+".pil.z.%f.%s" % (args.z, args.format))
        

if __name__ == "__main__":
    main()