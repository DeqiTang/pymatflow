#!/usr/bin/env python

import os
import sys
import numpy as np
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
    elif filepath.split(".")[-1] == "lammps" or filepath.split(".")[-1] == "lmp":
        from pymatflow.structure.crystal import crystal
        import pymatflow.third.aseio as aseio
        a = crystal()
        a.cell, a.atoms = aseio.read_lammps_data(filepath)
    else:
        pass
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
    elif filepath.split(".")[-1] == "lammps" or filepath.split(".")[-1] == "lmp":
        import pymatflow.third.aseio as aseio
        aseio.write_lammps_data(cell=structure.cell, atoms=structure.atoms, filepath=filepath)            
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

    subparser.add_argument("--fix", help="list of fixed atoms, index start from 1, have privilege over --around-z", nargs='+', type=int, default=None)

    subparser.add_argument("--around-z", type=float, nargs=3, default=None,
        help="select atoms around specified z in Angstrom with tolerance, like this --around-z 10 -0.5 0.5")

    subparser.add_argument("--color", type=str, default="white",
        choices=["red", "green", "blue", "white"],
        help="select color to color the fix atoms in xsd file, can be: red green blue and white")

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

    subparser.add_argument("--thick", type=float, default=10,
            help="thickness of the vacuum layer, in unit of Angstrom, default is 10")

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
            
    subparser.add_argument("--precision", type=float, default=1.0e-8,
            help="a value that is less than 1 and infinitely close to 1 used to judge whether one atom is in another periodic of the redefined cell")
            
    # ---------------------------------------------------------------------------------
    # cleave surface
    # ---------------------------------------------------------------------------------
    subparser = subparsers.add_parser("cleave", help="cleave surface")

    subparser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    subparser.add_argument("-o", "--output", type=str, required=True,
            help="output structure file")

    subparser.add_argument("--direction", type=int, nargs=3, default=[0, 0, 1],
            help="direction of the surface plane to cleave")            
            
    subparser.add_argument("--thick", type=float,
            help="thickness of the vacuum layer, in unit of Angstrom, default is 10")
            
    subparser.add_argument("--precision", type=float, default=1.0e-8,
            help="a value that is large than 0 and infinitely close to 0 used to judge whether one atom is in another periodic of the redefined cell used in cleave surface")
            
    # ---------------------------------------------------------------------------------
    # merge layers | ab plane
    # ---------------------------------------------------------------------------------
    subparser = subparsers.add_parser("merge", help="merge layers | ab plane")

    subparser.add_argument("-i", "--input", type=str, nargs=2, required=True,
            help="input structure files")

    subparser.add_argument("-o", "--output", type=str, required=True,
            help="output structure file")

    #subparser.add_argument("--direction", type=int, nargs=3, default=[0, 0, 1],
    #        help="direction of the surface plane to cleave")            
            
    subparser.add_argument("--usecell", type=str, default="average",
            choices=["1", "2", "average"],
            help="use cell of structure 1 or 2 , otherwise average by default")            
            
    subparser.add_argument("--thick", type=float,
            help="thickness of the vacuum layer, in unit of Angstrom, default is 10")
            
    subparser.add_argument("--distance", type=float,
            help="distance between the layer, in unit of Angstrom, default is 3.4")
            
    # ---------------------------------------------------------------------------------
    # nanotube builder
    # ---------------------------------------------------------------------------------
    subparser = subparsers.add_parser("tube", help="nanotube along b direction(a must be perpendicular to b and ab is the surface plane)")

    subparser.add_argument("-i", "--input", type=str, required=True,
            help="input structure files")

    subparser.add_argument("-o", "--output", type=str, required=True,
            help="output structure file")
            
    subparser.add_argument("--plane", type=int, default=1,
            help="on which plane to add vacuum layer. 1: ab, 2: ac, 3: bc")
            
    subparser.add_argument("--axis", type=str, default="b",
            choices=["a", "b", "c"],
            help="build nanotube along an axis parallel to axis specified")
     
    # -----------------------------------------------------------------------------------
    # set frac within zero and one
    # ------------------------------------------------------------------------------------
    subparser = subparsers.add_parser("std", help="set fractional coordinates within zero and one")
    
    subparser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")
            
    subparser.add_argument("-o", "--output", type=str, required=True,
            help="output structure file")
     
    # ------------------------------------------------------------------------------------
    # generate series of cell volume changed structures
    # ------------------------------------------------------------------------------------
    subparser = subparsers.add_parser("cv", help="generate series of cell volume changed structures")

    subparser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    subparser.add_argument("-d", "--directory", type=str, default="./",
            help="directory to put the generated structures")

    subparser.add_argument("--range", type=float, nargs=3, default=[0.95, 1.05, 0.01],
            help="cell volume change ratio, default is [0.95, 1.05, 0.01]")

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
        a = read_structure(filepath=args.input)        
        if args.fix != None:
            fix = args.fix
        elif args.around_z != None:
            atoms_index_from_1 = []
            for i in range(len(a.atoms)):
                if a.atoms[i].z > (args.around_z[0] + args.around_z[1]) and a.atoms[i].z < (args.around_z[0] + args.around_z[2]):
                    atoms_index_from_1.append(i+1)
            fix = atoms_index_from_1
        else:
            fix = []
            
        if args.output.split(".")[-1] == "xyz":
            fix_str = ""
            for i in fix:
                fix_str += "%d " % i
            os.system("xyz-fix-atoms.py -i %s -o %s --fix %s" % (args.input, args.output, fix_str))
        elif os.path.basename(args.output) == "POSCAR":
            from pymatflow.vasp.base.poscar import vasp_poscar
            for i in fix:
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
        # output an xsd file with fixed atoms colored specifically so that user can check the atoms fixed
        from xml.etree.ElementTree import parse
        os.system("mkdir -p /tmp/structflow/fix")
        write_structure(a, filepath="/tmp/structflow/fix/tmp.xsd")
        # read xsd file
        xsd = parse("/tmp/structflow/fix/tmp.xsd")
    
        # ID of Atom3D in xsd file start from 4
        imap = xsd.getroot().find("AtomisticTreeRoot").find("SymmetrySystem").find("MappingSet").find("MappingFamily").find("IdentityMapping")
        atoms = imap.findall("Atom3d")
        if args.color == "white":
            RGB = [255, 255, 255]
        elif args.color == "red":
            RGB = [255, 0, 0]
        elif args.color == "green":
            RGB = [0, 255, 0]
        elif args.color == "blue":
            RGB = [0, 0, 255]
        else:
            RGB = [255, 255, 255] # default
            
        for i in fix:
            atoms[i-1].set("Color", "%f, %f, %f, %f" % (RGB[0], RGB[1], RGB[2], 1))

        # write xsd file
        xsd.write(args.input+".coloring.atoms.fixed.xsd")    
        

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

        vacuum_layer(a, plane=args.plane, thickness=args.thick if args.thick != None else 10.0)
        
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

        redefined = redefine_lattice(structure=a, a=args.a, b=args.b, c=args.c, precision=args.precision)
        
        # output structure
        write_structure(structure=redefined, filepath=args.output)        
    elif args.driver == "cleave":
        from pymatflow.structure.tools import cleave_surface
        a = read_structure(filepath=args.input) 
        print("=======================================================================\n")
        print("                       structflow\n")
        print("-----------------------------------------------------------------------\n")
        print("you are trying to cleave the surface of (%d, %d, %d)\n" % (args.direction[0], args.direction[1], args.direction[2]))            
        print("from %s\n" % args.input)
        print("\n")
        print("the output structure file is -> %s\n" % args.output)

        cleaved = cleave_surface(structure=a, direction=args.direction, thickness=args.thick if args.thick != None else 10.0, precision=args.precision)
        
        # output structure
        write_structure(structure=cleaved, filepath=args.output)           
    elif args.driver == "merge":
        from pymatflow.structure.tools import merge_layers
        a_list = []
        for i in range(2):
            a_list.append(read_structure(filepath=args.input[i]))
        print("=======================================================================\n")
        print("                       structflow\n")
        print("-----------------------------------------------------------------------\n")
        print("you are trying to merge layers on ab plane\n")            
        print("from %s\n" % (args.input[0]))
        print("and %s\n" % (args.input[1]))
        print("\n")
        print("the output structure file is -> %s\n" % args.output)

        if args.usecell == "1":
            usecell = 1
        elif args.usecell == "2":
            usecell = 2
        else:
            usecell = "average"
        merged = merge_layers(structure1=a_list[0], structure2=a_list[1], use_cell=usecell, distance=args.distance if args.distance != None else 3.4, thickness=args.thick if args.thick != None else 10.0)
        
        # output structure
        write_structure(structure=merged, filepath=args.output)                 
    elif args.driver == "tube":
        a = read_structure(filepath=args.input)
        if args.plane == 1:
            plane = "ab"
        elif args.plane == 2:
            plane = "ac"
        elif args.plane == 3:
            plane = "bc"        
        print("=======================================================================\n")
        print("                       structflow\n")
        print("-----------------------------------------------------------------------\n")
        print("you are trying to build nanotube of %s plane along %s vector\n" % (plane, args.axis))            
        print("from %s\n" % (args.input))
        print("the output structure file is -> %s\n" % args.output)
        tube = None
        if plane == "ab":
            from pymatflow.structure.tools import build_nanotube_ab        
            if args.axis not in "ab":
                print("building nanotube of ab plane along axis parallel to c is unphysical!!!\n")
                sys.exit()
            else:
                tube = build_nanotube_ab(structure=a, axis=args.axis)
        if plane == "ac":
            from pymatflow.structure.tools import build_nanotube_ac
            if args.axis not in "ac":
                print("building nanotube of ac plane along axis parallel to b is unphysical!!!\n")
                sys.exit()                
            else:
                tube = build_nanotube_ac(structure=a, axis=args.axis)
        if plane == "bc":
            from pymatflow.structure.tools import build_nanotube_bc
            if args.axis not in "bc":
                print("building nanotube of bc plane along axis parallel to a is unphysical!!!\n")
                sys.exit()                
            else:
                tube = build_nanotube_bc(structure=a, axis=args.axis)

        # output structure
        if tube != None:
            write_structure(structure=tube, filepath=args.output)             
    elif args.driver == "std":
        from pymatflow.structure.tools import set_frac_within_zero_and_one
        a = read_structure(filepath=args.input)
        print("=======================================================================\n")
        print("                       structflow\n")
        print("-----------------------------------------------------------------------\n")
        print("you are trying to set fractional coords within 0 and 1\n")            
        print("from %s\n" % (args.input))
        print("\n")
        print("the output structure file is -> %s\n" % args.output)
        
        normalized = set_frac_within_zero_and_one(structure=a)
        # output structure
        write_structure(structure=normalized, filepath=args.output)           
    elif args.driver == "cv":
        from pymatflow.structure.crystal import crystal
        from pymatflow.base.atom import Atom
        a = read_structure(filepath=args.input)
        print("=======================================================================\n")
        print("                       structflow\n")
        print("-----------------------------------------------------------------------\n")
        print("you are trying to get a series of structure with different volume\n")            
        print("from %s\n" % (args.input))
        print("\n")
        print("the output dir for structure file is -> %s\n" % args.directory)

        # now calc the fractional coordinates
        atoms_frac = []
        latcell = np.array(a.cell)
        convmat = np.linalg.inv(latcell.T)
        for i in range(len(a.atoms)):
            atom = []
            atom.append(a.atoms[i].name)
            atom = atom + list(convmat.dot(np.array([a.atoms[i].x, a.atoms[i].y, a.atoms[i].z])))
            atoms_frac.append(atom)
        #
        out = crystal()
        os.system("mkdir -p %s" % args.directory)
        for i, ratio_v in enumerate(np.arange(args.range[0], args.range[1], args.range[2])):
            ratio = np.power(ratio_v, 1/3)

            # now convert coord of atom in atoms_frac_within_new_cell to cartesian
            out.atoms = []
            out.cell = (np.array(a.cell) * ratio).tolist()
            latcell = np.array(out.cell)
            convmat_frac_to_cartesian = latcell.T
            for atom in atoms_frac:
                cartesian = list(convmat_frac_to_cartesian.dot(np.array([atom[1], atom[2], atom[3]])))
                out.atoms.append(Atom(name=atom[0], x=cartesian[0], y=cartesian[1], z=cartesian[2]))
            output_name = ".".join(os.path.basename(args.input).split(".")[:-1] + ["%d" % i, "cif"])
            write_structure(out, filepath=os.path.join(args.directory, output_name))
            #            
        with open(os.path.join(args.directory, "log.txt"), 'w') as fout:
            fout.write("# index\tratio_v\tvolume(Angstrom^3)\n")
            for i, ratio_v in enumerate(np.arange(args.range[0], args.range[1], args.range[2])):
                ratio = np.power(ratio_v, 1/3)
                cell_now = (np.array(a.cell) * ratio).tolist()
                fout.write("%d\t%f\t%f\n" % (i, ratio_v, np.linalg.det(cell_now)))
    # --------------------------------------------------------------------------



if __name__ == "__main__":
    main()
