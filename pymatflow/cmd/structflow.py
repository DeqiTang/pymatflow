#!/usr/bin/env python

import os
import sys
import argparse


def read_structure(filepath):
    """
    read in input structure

    :param filepath: file path for the input structure file
        it will judge the file type by the suffix of the file

    :return a: an instance of pymatflow.structure.crystal
    """
    if filepath.split(".")[-1] == "xyz":
        from pymatflow.structure.crystal import crystal
        a = crystal()
        a.from_xyz_file(filepath)
    elif filepath.split(".")[-1] == "cif":
        from pymatflow.structure.crystal import crystal
        import pymatflow.third.aseio as aseio
        a = crystal()
        a.cell, a.atoms = aseio.read_cif(filepath)
    elif filepath.split(".")[-1] == "xsd":
        from pymatflow.structure.crystal import crystal
        import pymatflow.third.aseio as aseio
        a = crystal()
        a.cell, a.atoms = aseio.read_xsd(filepath)
    elif filepath.split(".")[-1] == "xsf":
        from pymatflow.structure.crystal import crystal
        import pymatflow.third.aseio as aseio
        a = crystal()
        a.cell, a.atoms = aseio.read_xsf(filepath)
    elif os.path.basename(filepath) == "POSCAR" or os.path.basename(filepath) == "CONTCAR":
        from pymatflow.structure.crystal import crystal
        import pymatflow.third.aseio as aseio
        a = crystal()
        a.cell, a.atoms = aseio.read_poscar(filepath)
    return a


def write_structure(structure, filepath, frac=1):
    """
    write structure to file
    :param structure: an instance of pymatflow.structure.crystal

    :param filepath: file path for the output structure file
        it will judge the file type by the suffix
    :param frac: output fractional coordinates, currently only used by POSCAR/CONTCAR
        1(default): use faractional, 0: use cartesian
    """
    if filepath.split(".")[-1] == "xyz":
        structure.write_xyz(filepath=filepath)
    elif filepath.split(".")[-1] == "cif":
        import pymatflow.third.aseio as aseio
        aseio.write_cif(cell=structure.cell, atoms=structure.atoms, filepath=filepath)
    elif filepath.split(".")[-1] == "xsd":
        import pymatflow.third.aseio as aseio
        aseio.write_xsd(cell=structure.cell, atoms=structure.atoms, filepath=filepath)
    elif filepath.split(".")[-1] == "xsf":
        import pymatflow.third.aseio as aseio
        aseio.write_xsf(cell=structure.cell, atoms=structure.atoms, filepath=filepath)
    elif filepath.split(".")[-1] == "cube":
        import pymatflow.third.aseio as aseio
        aseio.write_cube(cell=structure.cell, atoms=structure.atoms, filepath=filepath)
    elif os.path.basename(filepath) == "POSCAR" or os.path.basename(filepath) == "CONTCAR":
        from pymatflow.structure.crystal import crystal
        from pymatflow.vasp.base.poscar import vasp_poscar
        #import pymatflow.third.aseio as aseio
        poscar = vasp_poscar()
        poscar.xyz.cell = structure.cell
        poscar.xyz.atoms = structure.atoms
        poscar.xyz.natom = len(poscar.xyz.atoms)
        poscar.xyz.set_species_number() # needed for poscar output
        with open(filepath, 'w') as fout:
            #poscar.to_poscar(fout=fout, coordtype="Cartesian" if frac == 0 else "Direct")            
            if frac == 0:
                coordtype = "Cartesian"
            elif frac == 1:
                coordtype = "Direct"
            else:
                pass
            poscar.to_poscar(fout=fout, coordtype=coordtype)
    else:
        pass

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="driver", title="subcommands", description="choose one and only one subcommand")

    
    # --------------------------------------------------------------------------
    # supercell builder
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("supercell", help="using supercell subcommand")

    subparser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    subparser.add_argument("-o", "--output", type=str, required=True,
            help="output structure file")

    subparser.add_argument("-n", "--supern", nargs="+", type=int,
            help="bulid supern:[int, int, int] supercell")


    # --------------------------------------------------------------------------
    # fix atoms
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("fix", help="using fix subcommand")

    subparser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    subparser.add_argument("-o", "--output", type=str, required=True,
            help="output structure file")

    subparser.add_argument("--fix", help="list of fixed atoms", nargs='+', type=int)


    # --------------------------------------------------------------------------
    # convert file type
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("convert", help="using convert subcommand")

    subparser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    subparser.add_argument("-o", "--output", type=str, required=True,
            help="output structure file")

    # --------------------------------------------------------------------------
    # kpath
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("kpath", help="using kpath subcommand")

    subparser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    subparser.add_argument("--engine", type=str, default="seekpath",
            choices=["seekpath"],
            help="choose tool to generate kpath")

    subparser.add_argument("--kpath-file", type=str, default="kpath-from-seekpath.txt",
            help="the output kpoints file")


    # ---------------------------------------------------------------------------------
    # move atoms along one direction
    # ---------------------------------------------------------------------------------
    subparser = subparsers.add_parser("move", help="move atoms along one direction")

    subparser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    subparser.add_argument("-o", "--output", type=str, required=True,
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

    subparser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    subparser.add_argument("-o", "--output", type=str, required=True,
            help="output structure file")

    subparser.add_argument("--atoms", type=int, nargs="+",
            help="atoms to remove, index start from 1")

    subparser.add_argument("--elements", type=str, nargs="+",
            help="elements to remove")

    # ---------------------------------------------------------------------------------
    # vacuum layer
    # ---------------------------------------------------------------------------------
    subparser = subparsers.add_parser("vacuum", help="add vacuum layer")

    subparser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    subparser.add_argument("-o", "--output", type=str, required=True,
            help="output structure file")

    subparser.add_argument("--plane", type=int, default=1,
            help="on which plane to add vacuum layer. 1: ab, 2: ac, 3: bc")

    subparser.add_argument("--thick", type=float,
            help="thickness of the vacuum layer, in unit of Angstrom")

    # ---------------------------------------------------------------------------------
    # inverse atoms against geometric center
    # ---------------------------------------------------------------------------------
    subparser = subparsers.add_parser("inverse", help="inverse against geo center")

    subparser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    subparser.add_argument("-o", "--output", type=str, required=True,
            help="output structure file")

    subparser.add_argument("-c", "--center", type=str, default="cell",
            choices=["geo", "cell"],
            help="inversion center, can geo or cell")
            
    # ---------------------------------------------------------------------------------
    # redefine lattice
    # ---------------------------------------------------------------------------------
    subparser = subparsers.add_parser("redefine", help="redefine lattice")

    subparser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    subparser.add_argument("-o", "--output", type=str, required=True,
            help="output structure file")

    subparser.add_argument("-a", type=int, nargs=3, default=[1, 0, 0],
            help="a from old a b c")
            
    subparser.add_argument("-b", type=int, nargs=3, default=[0, 1, 0],
            help="b from old a b c")            
            
    subparser.add_argument("-c", type=int, nargs=3, default=[0, 0, 1],
            help="c from old a b c")            
            
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
        
        a = read_structure(filepath=args.input)
        supercell = a.build_supercell(args.supern)
        new_structure = crystal()
        new_structure.get_cell_atoms(cell=supercell["cell"], atoms=supercell["atoms"])
        write_structure(structure=new_structure, filepath=args.output)
        
        print("=========================================================\n")
        print("              structflow supercell builder\n")
        print("---------------------------------------------------------\n")
        print("you are trying to bulid supercell from %s\n" % args.input)
        print("the output structure file is -> %s\n" % args.output)

    elif args.driver == "fix":
        # can only write xyz and poscar file
        
        if args.output.split(".")[-1] == "xyz":
            fix_str = ""
            for i in args.fix:
                fix_str += "%d " % i
            os.system("xyz-fix-atoms.py -i %s -o %s --fix %s" % (args.input, args.output, fix_str))
        elif os.path.basename(args.output) == "POSCAR":
            from pymatflow.vasp.base.poscar import vasp_poscar
            a = read_structure(filepath=args.input)
            for i in args.fix:
                a.atoms[i-1].fix = [True, True, True]
            poscar = vasp_poscar()
            poscar.xyz.cell = a.cell
            poscar.xyz.atoms = a.atoms
            poscar.xyz.natom = len(poscar.xyz.atoms)
            poscar.xyz.set_species_number() # needed for poscar output
            poscar.selective_dynamics = True
            with open(args.output, 'w') as fout:
                poscar.to_poscar(fout=fout, coordtype="Direct")
        else:
            print("===============================================================\n")
            print("                      WARNING !!!\n")
            print("---------------------------------------------------------------\n")
            print("structflow fix now only supports write of xyz and POSCAR\n")
            sys.exit(1)        

    elif args.driver == "convert":
        # will convert file type according to the suffix of the specified input and output file

        a = read_structure(filepath=args.input)
        write_structure(structure=a, filepath=args.output)


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
        a = read_structure(filepath=args.input)
        # move atoms
        print("=========================================================\n")
        print("                   structflow\n")
        print("----------------------------------------------------------\n")
        print("you are trying to move atoms:\n")
        print(args.atoms)
        for i in args.atoms:
            print("%d -> %s\n" % (i, a.atoms[i-1].name))
        print("\n")
        print("along direction:\n")
        print(args.direction)
        print("\n")
        print("by length of -> %f, in unit of Angstrom\n" % args.disp)
        move_along(a, atoms_to_move=[i-1 for i in args.atoms], direc=args.direction, disp=args.disp)
        
        # output structure
        write_structure(structure=a, filepath=args.output)

    elif args.driver == "remove":
        from pymatflow.structure.tools import remove_atoms
        a = read_structure(filepath=args.input)
        # remove atoms
        print("=======================================================================\n")
        print("                       structflow\n")
        print("-----------------------------------------------------------------------\n")
        print("you are trying to remove from %s the following list of atoms:\n" % args.input)
        print(args.atoms)
        if args.atoms != None:
            for i in args.atoms:
                print("%d -> %s\n" % (i, a.atoms[i-1].name))
        else:
            pass
        print("\n")
        print("also the following elements will be removed:\n")
        print(args.elements)
        print("the output structure file is -> %s\n" % args.output)

        if args.atoms != None:
            remove_atoms(a, atoms_to_remove=[i-1 for i in args.atoms])
       
        # we should first remove atoms specified by args.atoms
        # and remove atoms specified by args.elements
        # as remove atom will change the index of atom
        if args.elements != None:
            element_atoms_to_remove = []
            for i in range(len(a.atoms)):
                if a.atoms[i].name in args.elements:
                    element_atoms_to_remove.append(i)
            remove_atoms(a, atoms_to_remove=element_atoms_to_remove)

        # output structure
        write_structure(structure=a, filepath=args.output)

    elif args.driver == "vacuum":
        from pymatflow.structure.tools import vacuum_layer
        a = read_structure(filepath=args.input) 
        # add vacuum layer
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
        write_structure(structure=a, filepath=args.output)
    elif args.driver == "inverse":
        from pymatflow.structure.tools import inverse_geo_center
        from pymatflow.structure.tools import inverse_cell_center
        a = read_structure(filepath=args.input) 
        print("=======================================================================\n")
        print("                       structflow\n")
        print("-----------------------------------------------------------------------\n")
        if args.center == "geo":
            print("you are trying to inverse the system against the geometric center\n")
        elif args.center == "cell":
            print("you are trying to inverse the system against the cell center\n")            
        print("from %s\n" % args.input)
        print("\n")
        print("the output structure file is -> %s\n" % args.output)

        if args.center == "geo":
            inverse_geo_center(a)
        elif args.center == "cell":
            inverse_cell_center(a)
        
        # output structure
        write_structure(structure=a, filepath=args.output)      
    elif args.driver == "redefine":
        from pymatflow.structure.tools import redefine_lattice
        a = read_structure(filepath=args.input) 
        print("=======================================================================\n")
        print("                       structflow\n")
        print("-----------------------------------------------------------------------\n")
        print("you are trying to redefine the lattice\n")            
        print("from %s\n" % args.input)
        print("\n")
        print("the output structure file is -> %s\n" % args.output)

        redefined = redefine_lattice(structure=a, a=args.a, b=args.b, c=args.c)
        
        # output structure
        write_structure(structure=redefined, filepath=args.output)              
    # --------------------------------------------------------------------------



if __name__ == "__main__":
    main()
