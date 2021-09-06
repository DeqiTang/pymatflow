#!/usr/bin/env python

import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

from pymatflow.cmd.structflow import read_structure
from pymatflow.cmd.structflow import write_structure

class VaspCHG:
    def __init__(self, filepath=None):
        self.structure = None
        if filepath != None:
            self.get_chg(filepath)

    def get_chg(self, filepath):
        """
        :param filepath: the *CHG* file path
        """
        with open(filepath, "r") as fin:
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
                    
        os.system("mkdir -p /tmp/pymatflow/charge/chg_vasp")
        with open("/tmp/pymatflow/charge/chg_vasp/POSCAR", "w") as fout:
            for i in range(first_blank_line):
                fout.write(chg[i])
                
        self.structure = read_structure("/tmp/pymatflow/charge/chg_vasp/POSCAR")
        
        is_orthogonal = False
        cos_ab = np.dot(np.array(self.structure.cell[0]), np.array(self.structure.cell[1]))
        cos_ac = np.dot(np.array(self.structure.cell[0]), np.array(self.structure.cell[2]))
        cos_bc = np.dot(np.array(self.structure.cell[1]), np.array(self.structure.cell[2]))
        if cos_ab == cos_ac == cos_bc == 0:
            is_orthogonal = True
        self.is_orthogonal = is_orthogonal

        self.ngxf = int(chg[first_blank_line+1].split()[0])
        self.ngyf = int(chg[first_blank_line+1].split()[1])
        self.ngzf = int(chg[first_blank_line+1].split()[2])    
        
        if first_augmentation_line == None:
            #data = np.loadtxt(chg[first_blank_line+2:])
            tmp_str = "".join(chg[first_blank_line+2:])
            data = np.fromstring(tmp_str, sep="\n")
        else:
            #data = np.loadtxt(chg[first_blank_line+2:first_augmentation_line])
            tmp_str = "".join(chg[first_blank_line+2:first_augmentation_line])
            data = np.fromstring(tmp_str, sep="\n")
            
        self.data = data.reshape(self.ngzf, self.ngyf, self.ngxf)
        # charge data in cube file is in shape (ngridx, ngridy, ngridz)
        # while charge in *CHG* file is in shape (ngzf, ngyf, ngxf)
        # they are different!    

        # the unit of value is actually not physical now!
        self.cell_volume = np.dot(np.cross(np.array(self.structure.cell[0]), np.array(self.structure.cell[1])), np.array(self.structure.cell[2]))
        self.cell_volume_per_unit = self.cell_volume / (self.ngzf * self.ngyf * self.ngxf)
        # value in Vasp *CHG* are \rho(r)_of_electrons * Volume_of_cell, so we should divide it by cell_volume here and time it with cell_volume_per_unit
        # to get the number of electrons per divided unit
        self.total_electrons = np.sum(self.data) / self.cell_volume * self.cell_volume_per_unit
        #
        
        #print("======================================================\n")
        #print("           Information collected\n")
        #print("------------------------------------------------------\n")
        #print("cell volume: %f (A^3)\n" % self.cell_volume)
        #print("total electrons: %f\n" % self.total_electrons)

    
    def plot_grayscale_z(self, z=1, output_prefix="chg"):
        """
        # -------------------------------------------------------
        # gray scale image only for z direction
        # may not work for triclinic and monoclinic crystal system
        # -------------------------------------------------------
        :param z: a value between 0 and 1, indicating
                height of in z direction to print the plot
        :param output_prefix: prefix of the output image file name
        """
        zi = int((self.data.shape[0]-1) * z)
        #img = self.data[i, ::-1, ::]
        img = self.data[zi, ::, ::]
        img = (img-img.min()) / (img.max() - img.min()) * 255
        # need to do a transform when the cell is not Orthorhombic
        # skew the image
        a = np.array(self.structure.cell[0])
        b = np.array(self.structure.cell[1])
        cosangle = a.dot(b)/(np.linalg.norm(a) * np.linalg.norm(b))
        angle = np.arccos(cosangle) * 180 / np.pi        
        ax = plt.axes() #plt.figure()
        n1 = self.ngxf
        n2 = self.ngyf
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
        plt.savefig(output_prefix+"-z-%f.grayscale.plot.%s" % (z, "png"))
        #plt.imsave(args.output+"-z-%f.grayscale.image.png" % args.z, arr=img, cmap="gray") # img is not Affine transoformed here!!!
        plt.close()
        
    def plot_contour_2d(self, z=1, levels=10, cmap="gray", output_prefix="chg"):
        """
        # -----------------------------------------------------------------------------
        # 2D contour plot
        #------------------------------------------------------------------------------
        :param z: a value between 0 and 1, indicating
            height of in z direction to print the plot
        :param levels: levels of the color map or color bar
        :param cmap: color map style
        :param output_prefix: prefix of the output image file name
        """
        zi = int((self.data.shape[0]-1) * z)        
        nx = np.linspace(0, 1, self.ngxf)
        ny = np.linspace(0, 1, self.ngyf)
        X, Y = np.meshgrid(nx, ny) # now this Mesh grid cannot be used directly, we have to calc the real x y for it
        for xi in range(len(nx)):
            for yi in range(len(ny)):
                X[yi, xi] = self.structure.cell[0][0] * nx[xi] + self.structure.cell[1][0] * ny[yi]
                Y[yi, xi] = self.structure.cell[0][1] * nx[xi] + self.structure.cell[1][1] * ny[yi]
    
        # export data
        density_z = self.data[zi, :, :] / self.cell_volume
        with open(output_prefix+'.slice-z.%f.data' % z, "w") as fout:
            fout.write("# x y val(e/Angstrom^3)\n")
            for xi in range(len(nx)):
                for yi in range(len(ny)):
                    fout.write("%f %f %f\n" % (X[xi, yi], Y[xi, yi], density_z[xi, yi]))

        Z = self.data[zi, :, :]
        Z = (Z-Z.min()) / (Z.max() - Z.min()) * 255
        # fill color, three color are divided into three layer(6)
        # cmap = plt.cm.hot means using thermostat plot(graduated red yellow)
        #cset = plt.contourf(X, Y, Z, levels=args.levels, cmap=plt.cm.hot)
        #cset = plt.contourf(X, Y, Z, levels=args.levels, cmap=plt.cm.gray)
        cset = plt.contourf(X, Y, Z, levels=levels, cmap=cmap)
        contour = plt.contour(X, Y, Z, levels=[20, 40], colors='k')
        plt.colorbar(cset)
        plt.autoscale()
        plt.tight_layout()
        plt.axis("equal") # set axis equally spaced
        #plt.show()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig(output_prefix+".2d-clice-z-%f.%s" % (z, "png"))
        plt.close()
        
    
    def plot_grayscale_orthogonal(self, z=1, output_prefix="chg"):
        """
        # -------------------------------------------------------
        # gray scale image only for orthogonal crystal system
        # -------------------------------------------------------
        :param z: a value between 0 and 1, indicating
                height of in z direction to print the plot
        :param output: prefix of the output image file name
        """
        zi = int((self.data.shape[0]-1) * z)
        #img = data[i, ::-1, ::]
        img = self.data[zi, ::, ::]
        img = (img-img.min()) / (img.max() - img.min()) * 255
        
        import PIL.Image as Image
        
        img = Image.fromarray(img.astype('uint8'))
        img = img.convert('L') # gryscale
        img.save(output_prefix+".pil.z.%f.%s" % (z, "png"))
        