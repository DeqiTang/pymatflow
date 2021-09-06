#!/usr/bin/env python

import argparse

from pymatflow.cmd.structflow import read_structure
from pymatflow.cmd.structflow import write_structure


def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str, required=True,
        help="input vasp *CHG* file")

    parser.add_argument("--output-structure", type=str, default="chg.cif",
        help="output stucture contained in PARCHG")
        
    parser.add_argument("-o", "--output", type=str, default="chg",
        help="prefix of the output image file name")         
    
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

    from pymatflow.charge.chg_vasp import VaspCHG

    vaspchg = VaspCHG()
    vaspchg.get_chg(chg_filepath)

    
    if 1 in args.output_option:
        # -------------------------------------------------------
        # gray scale image only for z direction
        # may not work for triclinic and monoclinic crystal system
        # -------------------------------------------------------
        vaspchg.plot_grayscale_z(z=args.z, output_prefix=args.output)
    
    if 2 in args.output_option:
        # -----------------------------------------------------------------------------
        # 2D contour plot
        #------------------------------------------------------------------------------
        vaspchg.plot_contour_2d(z=args.z, levels=args.levels, cmap=args.cmap, output_prefix=args.output)

    
    if (3 in args.output_option) and vaspchg.is_orthogonal == True:
        # -------------------------------------------------------
        # gray scale image only for orthogonal crystal system
        # -------------------------------------------------------
        vaspchg.plot_grayscale_orthogonal(z=args.z, output_prefix=args.output)

    write_structure(vaspchg.structure, filepath=args.output_structure)

if __name__ == "__main__":
    main()