#!/usr/bin/env python

import os
import sys
import argparse



def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="driver", title="subcommands", description="choose one and only one subcommand")

    # --------------------------------------------------------------------------
    # supercell builder
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("supercell", help="using supercell subcommand")

    subparser.add_argument("-i", "--input", type=str,
            help="input xyz file")

    subparser.add_argument("-o", "--output", type=str,
            help="output xyz file")

    subparser.add_argument("-n", "--supern", nargs="+", type=int,
            help="bulid supern:[int, int, int] supercell")


    # --------------------------------------------------------------------------
    # fix atoms
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("fix", help="using fix subcommand")

    subparser.add_argument("-i", "--input", help="intput xyz file", type=str)
    subparser.add_argument("-o", "--output", help="output xyz file", type=str)
    subparser.add_argument("--fix", help="list of fixed atoms", nargs='+', type=int)


    # --------------------------------------------------------------------------
    # convert file type
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("convert", help="using convert subcommand")

    subparser.add_argument("-i", "--input", type=str,
            help="input structure file")

    subparser.add_argument("-o", "--output", type=str,
            help="output structure file")

    # --------------------------------------------------------------------------
    # kpath
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("kpath", help="using kpath subcommand")

    subparser.add_argument("--engine", type=str, default="seekpath",
            choices=["seekpath"],
            help="choose tool to generate kpath")

    subparser.add_argument("-i", "--input", type=str, default=None,
            help="the input xyz structure file")

    subparser.add_argument("-o", "--output", type=str, default="kpath-from-seekpath.txt",
            help="the output kpoints file")


    # ---------------------------------------------------------------------------------
    # move atoms along one direction
    # ---------------------------------------------------------------------------------
    subparser = subparsers.add_parser("move", help="move atoms along one direction")

    subparser.add_argument("-i", "--input", type=str,
            help="input structure file")

    subparser.add_argument("-o", "--output", type=str,
            help="output structure file")

    subparser.add_argument("--atoms", type=int, nargs="+",
            help="atoms to move, index start from 1")

    subparser.add_argument("--direction", type=float, nargs=3,
            help="direction to move the atoms, in format of crystal orientation index")

    subparser.add_argument("--disp", type=float,
            help="displacement along the moving direction, in unit of Anstrom")

    # ---------------------------------------------------------------------------------
    # remove atoms
    # ---------------------------------------------------------------------------------
    subparser = subparsers.add_parser("remove", help="remove specified atoms")

    subparser.add_argument("-i", "--input", type=str,
            help="input structure file")

    subparser.add_argument("-o", "--output", type=str,
            help="output structure file")

    subparser.add_argument("--atoms", type=int, nargs="+",
            help="atoms to move, index start from 1")

    # ---------------------------------------------------------------------------------
    # vacuum layer
    # ---------------------------------------------------------------------------------
    subparser = subparsers.add_parser("vacuum", help="add vacuum layer")

    subparser.add_argument("-i", "--input", type=str,
            help="input structure file")

    subparser.add_argument("-o", "--output", type=str,
            help="output structure file")

    subparser.add_argument("--plane", type=int, default=1,
            help="on which plane to add vacuum layer. 1: ab, 2: ac, 3: bc")

    subparser.add_argument("--thick", type=float,
            help="thickness of the vacuum layer, in unit of Angstrom")

    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()

    # if no argument passed to matflow
    if len(sys.argv) == 1:
        # display help message when no args provided
        parser.print_help()
        sys.exit(1)



    if args.driver == "supercell":
        from pymatflow.base.xyz import base_xyz
        from pymatflow.structure.crystal import crystal
        xyz = base_xyz()
        xyz.get_xyz(args.input)
        structure = crystal()
        new_structure = crystal()
        structure.from_base_xyz(xyz)
        supercell = structure.build_supercell(args.supern)
        new_structure.get_cell_atoms(cell=supercell["cell"], atoms=supercell["atoms"])
        new_structure.to_base_xyz().to_xyz_file(args.output)
    elif args.driver == "fix":
        fix_str = ""
        for i in args.fix:
            fix_str += "%d " % i
        os.system("xyz-fix-atoms.py -i %s -o %s --fix %s" % (args.input, args.output, fix_str))
    elif args.driver == "convert":
        # will convert file type according to the suffix of the specified input and output file

        # input structure
        if args.input.split(".")[-1] == "xyz":
            from pymatflow.structure.crystal import crystal
            a = crystal()
            a.from_xyz_file(args.input)
        elif args.input.split(".")[-1] == "cif":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_cif(args.input)
        elif args.input.split(".")[-1] == "xsd":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_xsd(args.input)
        elif args.input.split(".")[-1] == "xsf":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_xsf(args.input)
        elif os.path.basename(args.input) == "POSCAR" or os.path.basename(args.input) == "CONTCAR":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_poscar(args.input)

        # output structure
        if args.output.split(".")[-1] == "xyz":
            a.write_xyz(filepath=args.output)
        elif args.output.split(".")[-1] == "cif":
            import pymatflow.third.aseio as aseio
            aseio.write_cif(cell=a.cell, atoms=a.atoms, filepath=args.output)
        elif args.output.split(".")[-1] == "xsd":
            import pymatflow.third.aseio as aseio
            aseio.write_xsd(cell=a.cell, atoms=a.atoms, filepath=args.output)
        elif args.output.split(".")[-1] == "xsf":
            import pymatflow.third.aseio as aseio
            aseio.write_xsf(cell=a.cell, atoms=a.atoms, filepath=args.output)
        elif args.output.split(".")[-1] == "cube":
            import pymatflow.third.aseio as aseio
            aseio.write_cube(cell=a.cell, atoms=a.atoms, filepath=args.output)
        else:
            pass

        print("=========================================================\n")
        print("              structflow convert\n")
        print("---------------------------------------------------------\n")
        print("with the help from ase.io\n")

    elif args.driver == "kpath":
        if args.engine == "seekpath":
            os.system("kpath-xyz-seekpath.py -i %s -o %s" % (args.input, args.output))
        else:
            pass
    elif args.driver == "move":
        from pymatflow.structure.tools import move_along
        # input structure
        if args.input.split(".")[-1] == "xyz":
            from pymatflow.structure.crystal import crystal
            a = crystal()
            a.from_xyz_file(args.input)
        elif args.input.split(".")[-1] == "cif":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_cif(args.input)
        elif args.input.split(".")[-1] == "xsd":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_xsd(args.input)
        elif args.input.split(".")[-1] == "xsf":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_xsf(args.input)
        elif os.path.basename(args.input) == "POSCAR" or os.path.basename(args.input) == "CONTCAR":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_poscar(args.input)
        # move atoms
        print("=========================================================\n")
        print("                   structflow\n")
        print("----------------------------------------------------------\n")
        print("you are trying to move atoms:\n")
        print(args.atoms)
        for i in args.atoms:
            print("%d -> %s\n" % (i, a.atoms[i].name))
        print("\n")
        print("along direction:\n")
        print(args.direction)
        print("\n")
        print("by length of -> %f, in unit of Angstrom\n" % args.disp)
        move_along(a, atoms_to_move=[i-1 for i in args.atoms], direc=args.direction, disp=args.disp)
        
        # output structure
        if args.output.split(".")[-1] == "xyz":
            a.write_xyz(filepath=args.output)
        elif args.output.split(".")[-1] == "cif":
            import pymatflow.third.aseio as aseio
            aseio.write_cif(cell=a.cell, atoms=a.atoms, filepath=args.output)
        elif args.output.split(".")[-1] == "xsd":
            import pymatflow.third.aseio as aseio
            aseio.write_xsd(cell=a.cell, atoms=a.atoms, filepath=args.output)
        elif args.output.split(".")[-1] == "xsf":
            import pymatflow.third.aseio as aseio
            aseio.write_xsf(cell=a.cell, atoms=a.atoms, filepath=args.output)
        elif args.output.split(".")[-1] == "cube":
            import pymatflow.third.aseio as aseio
            aseio.write_cube(cell=a.cell, atoms=a.atoms, filepath=args.output)
        else:
            pass        
    elif args.driver == "remove":
        from pymatflow.structure.tools import remove_atoms
        # input structure
        if args.input.split(".")[-1] == "xyz":
            from pymatflow.structure.crystal import crystal
            a = crystal()
            a.from_xyz_file(args.input)
        elif args.input.split(".")[-1] == "cif":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_cif(args.input)
        elif args.input.split(".")[-1] == "xsd":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_xsd(args.input)
        elif args.input.split(".")[-1] == "xsf":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_xsf(args.input)
        elif os.path.basename(args.input) == "POSCAR" or os.path.basename(args.input) == "CONTCAR":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_poscar(args.input)
        
        # remove atoms
        print("=======================================================================\n")
        print("                       structflow\n")
        print("-----------------------------------------------------------------------\n")
        print("you are trying to remove from %s the following list of atoms:\n" % args.input)
        print(args.atoms)
        for i in args.atoms:
            print("%d -> %s\n" % (i, a.atoms[i].name))
        print("\n")
        print("the output structure file is -> %s\n" % args.output)

        remove_atoms(a, atoms_to_remove=[i-1 for i in args.atoms])
        
        # output structure
        if args.output.split(".")[-1] == "xyz":
            a.write_xyz(filepath=args.output)
        elif args.output.split(".")[-1] == "cif":
            import pymatflow.third.aseio as aseio
            aseio.write_cif(cell=a.cell, atoms=a.atoms, filepath=args.output)
        elif args.output.split(".")[-1] == "xsd":
            import pymatflow.third.aseio as aseio
            aseio.write_xsd(cell=a.cell, atoms=a.atoms, filepath=args.output)
        elif args.output.split(".")[-1] == "xsf":
            import pymatflow.third.aseio as aseio
            aseio.write_xsf(cell=a.cell, atoms=a.atoms, filepath=args.output)
        elif args.output.split(".")[-1] == "cube":
            import pymatflow.third.aseio as aseio
            aseio.write_cube(cell=a.cell, atoms=a.atoms, filepath=args.output)
        else:
            pass      
    elif args.driver == "vacuum":
        from pymatflow.structure.tools import vacuum_layer
        # input structure
        if args.input.split(".")[-1] == "xyz":
            from pymatflow.structure.crystal import crystal
            a = crystal()
            a.from_xyz_file(args.input)
        elif args.input.split(".")[-1] == "cif":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_cif(args.input)
        elif args.input.split(".")[-1] == "xsd":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_xsd(args.input)
        elif args.input.split(".")[-1] == "xsf":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_xsf(args.input)
        elif os.path.basename(args.input) == "POSCAR" or os.path.basename(args.input) == "CONTCAR":
            from pymatflow.structure.crystal import crystal
            import pymatflow.third.aseio as aseio
            a = crystal()
            a.cell, a.atoms = aseio.read_poscar(args.input)
        
        # remove atoms
        print("=======================================================================\n")
        print("                       structflow\n")
        print("-----------------------------------------------------------------------\n")
        if args.plane == 1:
            plane = "ab"
        elif args.plane == 2:
            plane = "ac"
        elif args.plane == 3:
            plane = "bc"
        print("you are trying to add vacuum layer of %f Angstrom on %s plane\n" % (args.thick, plane))
        print("from %s\n" % args.input)
        print("\n")
        print("the output structure file is -> %s\n" % args.output)

        vacuum_layer(a, plane=args.plane, thickness=args.thick)
        
        # output structure
        if args.output.split(".")[-1] == "xyz":
            a.write_xyz(filepath=args.output)
        elif args.output.split(".")[-1] == "cif":
            import pymatflow.third.aseio as aseio
            aseio.write_cif(cell=a.cell, atoms=a.atoms, filepath=args.output)
        elif args.output.split(".")[-1] == "xsd":
            import pymatflow.third.aseio as aseio
            aseio.write_xsd(cell=a.cell, atoms=a.atoms, filepath=args.output)
        elif args.output.split(".")[-1] == "xsf":
            import pymatflow.third.aseio as aseio
            aseio.write_xsf(cell=a.cell, atoms=a.atoms, filepath=args.output)
        elif args.output.split(".")[-1] == "cube":
            import pymatflow.third.aseio as aseio
            aseio.write_cube(cell=a.cell, atoms=a.atoms, filepath=args.output)
        else:
            pass                        
    # --------------------------------------------------------------------------



if __name__ == "__main__":
    main()
