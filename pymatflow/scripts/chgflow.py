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



class vasp_chg:
    def __init__(self, chg_filepath):
        #self.chg_file_path = chg_filepath
        self.read_chg(chg_filepath)

    def read_chg(self, chg_filepath):
        with open(chg_filepath, "r") as fin:
            self.lines = fin.readlines()
        
        for j in range(len(self.lines)):
            if len(self.lines[j].split()) == 0:
                first_blank_line = j
                break
        
        first_augmentation_line = None
        for j in range(len(self.lines)):
            if "augmentation" in self.lines[j]:
                first_augmentation_line = j
                break
                    
        os.system("mkdir -p /tmp/pymatflow/")
        with open("/tmp/pymatflow/POSCAR", "w") as fout:
            for j in range(first_blank_line):
                fout.write(self.lines[j])
                
        self.structure = read_structure("/tmp/pymatflow/POSCAR")
        
        # assume three *CHG* have the same ngxf and ngyf ngzf
        self.ngxf = int(self.lines[first_blank_line+1].split()[0])
        self.ngyf = int(self.lines[first_blank_line+1].split()[1])
        self.ngzf = int(self.lines[first_blank_line+1].split()[2])    
        

        if first_augmentation_line == None:
            tmp_str = "".join(self.lines[first_blank_line+2:])
            self.data = np.fromstring(tmp_str, sep="\n").reshape(self.ngzf, self.ngyf, self.ngxf)

        # the unit of value is actually not physical now!
        self.cell_volume = np.dot(np.cross(np.array(self.structure.cell[0]), np.array(self.structure.cell[1])), np.array(self.structure.cell[2]))
        self.cell_volume_per_unit = self.cell_volume / (self.ngzf * self.ngyf * self.ngxf)
        # value in Vasp *CHG* are \rho(r)_of_electrons * Volume_of_cell, so we should divide it by cell_volume here and time it with cell_volume_per_unit
        # to get the number of electrons per divided unit
        # total_electrons = np.sum(data) / cell_volume * cell_volume_per_unit
        #
        
        #print("======================================================\n")
        #print("           Information collected\n")
        #print("------------------------------------------------------\n")
        #print("cell volume: %f (A^3)\n" % cell_volume)
        #print("total electrons: %f\n" % total_electrons)


def main():
    parser = argparse.ArgumentParser()
    
    subparsers = parser.add_subparsers(dest="driver", title="subcommands", description="choose one and only one subcommand")
    
    subparser = subparsers.add_parser("1d", help="dimension reduction of CHG* to one dimensional")

    subparser.add_argument("-i", "--input", type=str, required=True,
        help="input vasp *CHG* file, -i *CHG* ")

    subparser.add_argument("--output-structure", type=str, default="chg",
        help="output stucture contained in *CHG*")

    subparser.add_argument("-o", "--output", type=str, default="chg",
        help="prefix of the output image file name")
    
    subparser.add_argument("--levels", type=int, default=10,
        help="levels of the color map or color bar")

    subparser.add_argument("-z", "--z", type=float, default=1,
        help="a value between 0 and 1, indicat height of in z direction to print the plot")
        
    subparser.add_argument("--cmap", type=str, default="gray",
        choices=["gray", "hot", "afmhot", "Spectral", "plasma", "magma", "hsv", "rainbow", "brg"])
        
    subparser.add_argument("--abscissa", type=str, nargs="+", default=["a", "b", "c"], 
        choices=["a", "b", "c"], 
        help="choose the direction to do the dimension reduction")

    #  add 
    subparser = subparsers.add_parser("add", help="add CHG*")

    subparser.add_argument("-i", "--input", type=str, nargs="+", required=True,
        help="input vasp *CHG* file, -i CHG1 CHG2 CHG3 ... ")

    subparser.add_argument("--output-structure", type=str, default="chg",
        help="output stucture contained in *CHG*")

    subparser.add_argument("-o", "--output", type=str, default="chg-merged",
        help="prefix of the output chg file name")

    # slice
    subparser = subparsers.add_parser("slice", help="slice CHG* to get 2d DATA")

    subparser.add_argument("-i", "--input", type=str, required=True,
        help="input vasp *CHG* file, -i *CHG* ")

    subparser.add_argument("--output-structure", type=str, default="chg-slice",
        help="output stucture contained in *CHG*")

    subparser.add_argument("-o", "--output", type=str, default="chg",
        help="prefix of the output image file name")
    
    subparser.add_argument("--levels", type=int, default=10,
        help="levels of the color map or color bar")

    subparser.add_argument("-z", "--z", type=float, default=1,
        help="a value between 0 and 1, indicat height of in z direction to print the plot")
        
    subparser.add_argument("--cmap", type=str, default="gray",
        choices=["gray", "hot", "afmhot", "Spectral", "plasma", "magma", "hsv", "rainbow", "brg"])
        
    subparser.add_argument("--abscissa", type=str, nargs="+", default=["a", "b", "c"], 
        choices=["a", "b", "c"], 
        help="choose the direction to do the dimension reduction")
    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()
    
    if args.driver == "1d":

        chg_filepath = args.input
        
        with open(chg_filepath, "r") as fin:
            chg = fin.readlines()
        
        for j in range(len(chg)):
            if len(chg[j].split()) == 0:
                first_blank_line = j
                break
        
        first_augmentation_line = None
        for j in range(len(chg)):
            if "augmentation" in chg[j]:
                first_augmentation_line = j
                break
                    
        os.system("mkdir -p /tmp/pymatflow/")
        with open("/tmp/pymatflow/POSCAR", "w") as fout:
            for j in range(first_blank_line):
                fout.write(chg[j])
                
        structure = read_structure("/tmp/pymatflow/POSCAR")
        write_structure(structure=structure, filepath=args.output_structure+".cif")

        a = np.linalg.norm(structure.cell[0])
        b = np.linalg.norm(structure.cell[1])
        c = np.linalg.norm(structure.cell[2])
        
        # assume three *CHG* have the same ngxf and ngyf ngzf
        ngxf = int(chg[first_blank_line+1].split()[0])
        ngyf = int(chg[first_blank_line+1].split()[1])
        ngzf = int(chg[first_blank_line+1].split()[2])    
        

        if first_augmentation_line == None:
            tmp_str = "".join(chg[first_blank_line+2:])
            data = np.fromstring(tmp_str, sep="\n").reshape(ngzf, ngyf, ngxf)
        #data = data.reshape(ngzf, ngyf, ngxf)

        
        # data dimension reduction
        # the unit of value is actually not physical now!
        cell_volume = np.dot(np.cross(np.array(structure.cell[0]), np.array(structure.cell[1])), np.array(structure.cell[2]))
        cell_volume_per_unit = cell_volume / (ngzf * ngyf * ngxf)
        
        # value in Vasp *CHG* are \rho(r)_of_electrons * Volume_of_cell, so we should divide it by cell_volume here and time it with cell_volume_per_unit
        # to get the number of electrons per divided unit
        total_electrons = np.sum(data) / cell_volume * cell_volume_per_unit
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
            for ci in range(data.shape[0]):
                tmp = 0
                for bi in range(data.shape[1]):
                    tmp += np.sum(data[ci, bi, :])
                nelect_ci = tmp * factor
                rho_line = nelect_ci / len_ci
                data_red_c.append(rho_line)
        if "b" in args.abscissa:
            factor = cell_volume_per_unit / cell_volume
            len_bi = b / ngyf    
            for bi in range(data.shape[1]):
                tmp = 0
                for ai in range(data.shape[2]):
                    tmp += np.sum(data[:, bi, ai])
                nelect_bi = tmp * factor
                rho_line = nelect_bi / len_bi                
                data_red_b.append(rho_line)
        if "a" in args.abscissa:
            factor = cell_volume_per_unit / cell_volume
            len_ai = a / ngxf        
            for ai in range(data.shape[2]):
                tmp = 0
                for ci in range(data.shape[0]):
                    tmp += np.sum(data[ci, :, ai])
                nelect_ai = tmp * factor
                rho_line = nelect_ai / len_ai                   
                data_red_a.append(rho_line)    

        # output the data and make the plot
        if "c" in args.abscissa:
            with open(args.output+".1d.c.data", 'w') as fout:
                fout.write("#c(angstrom) rho(e) (number of electron per Angstrom)\n")
                c_coord = np.linspace(0, c, len(data_red_c))
                for i in range(len(data_red_c)):
                    fout.write("%f %f\n" % (c_coord[i], data_red_c[i]))
            plt.plot(np.linspace(0, c, len(data_red_c)), data_red_c)                
            plt.ylabel(r"$\rho (e/\AA)$")
            plt.tight_layout()
            plt.savefig(args.output+".1d.c.png")
            plt.close()                
        if "b" in args.abscissa:
            with open(args.output+".1d.b.data", 'w') as fout:
                fout.write("#b(angstrom) rho(e) (number of electron per Angstrom)\n")
                b_coord = np.linspace(0, b, len(data_red_b))
                for i in range(len(data_red_b)):
                    fout.write("%f %f\n" % (b_coord[i], data_red_b[i]))        
            plt.plot(np.linspace(0, b, len(data_red_b)), data_red_b)    
            plt.ylabel(r"$\rho (e/\AA)$")              
            plt.tight_layout()
            plt.savefig(args.output+".1d.b.png")
            plt.close()                
        if "a" in args.abscissa:
            with open(args.output+".1d.a.data", 'w') as fout:
                fout.write("#a(angstrom) rho(e) (number of electron per Angstrom)\n")
                a_coord = np.linspace(0, a, len(data_red_a))
                for i in range(len(data_red_a)):
                    fout.write("%f %f\n" % (a_coord[i], data_red_a[i]))
            plt.plot(np.linspace(0, a, len(data_red_a)), data_red_a)                
            plt.ylabel(r"$\rho (e/\AA)$")            
            plt.tight_layout()
            plt.savefig(args.output+".1d.a.png")
            plt.close()

    elif args.driver == "add":
        chgs = []   
        for item in args.input:
            chgs.append(vasp_chg(item))
        pass
       
    elif args.driver == "slice":
        chg = vasp_chg(args.input)
        zi = int((chg.data.shape[0]-1) * args.z)
        img = chg.data[zi, ::, ::]
        img = (img-img.min()) / (img.max() - img.min())  * 255
        # need to do a transform when the cell is not Orthorhombic
        # skew the image
        a = np.array(chg.structure.cell[0])
        b = np.array(chg.structure.cell[1])
        cosangle = a.dot(b) / (np.linalg.norm(a) * np.linalg.norm(b))
        angle = np.arccos(cosangle) * 180 / np.pi
        ax = plt.axes()
        n1 = chg.ngxf
        n2 = chg.ngyf
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
        plt.savefig(args.output+"-z-%f.grayscale.png" % args.z)
        plt.close()
     
        # -----------------------------------------------------------------------------
        # 2D contour plot
        #------------------------------------------------------------------------------
    
        nx = np.linspace(0, 1, chg.ngxf)
        ny = np.linspace(0, 1, chg.ngyf)
        X, Y = np.meshgrid(nx, ny) # now this mesh grid cannot be used directly, we have to calc the real x y for it
        for xi in range(len(nx)):
            for yi in range(len(ny)):
                X[yi, xi] = chg.structure.cell[0][0] * nx[xi] + chg.structure.cell[1][0] * ny[yi]
                Y[yi, xi] = chg.structure.cell[0][1] * nx[xi] + chg.structure.cell[1][1] * ny[yi]
   
        density_z = chg.data[zi, :, :] / chg.cell_volume
        # export data
        with open(args.output+'.slice-z.%f.data' % args.z, "w") as fout:
            fout.write("# x y val(e/Angstrom^3)\n")
            for xi in range(len(nx)):
                for yi in range(len(ny)):
                    fout.write("%f %f %f\n" % (X[xi, yi], Y[xi, yi], density_z[xi, yi]))

        Z = chg.data[zi, :, :]
        Z = (Z-Z.min()) / (Z.max() - Z.min()) * 255
        # fill color, three color are divided into three layer(6)
        # cmap = plt.cm.hot means using thermostat plot(graduated red yellow)
        #cset = plt.contourf(x, y, z, levels=args.levels, cmap=plt.cm.hot)
        #cset = plt.contourf(x, y, z, levels=args.levels, cmap=plt.cm.gray)
        cset = plt.contourf(X, Y, Z, levels=args.levels, cmap=args.cmap)
        contour = plt.contour(X, Y, Z, levels=[20, 40], colors='k')
        plt.colorbar(cset)
        plt.autoscale()
        plt.tight_layout()
        plt.axis("equal") # set axis equally spaced
        #plt.show()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig(args.output+".2d-slice-z-%f.png" % args.z)
        plt.close()   

if __name__ == "__main__":
    main()
