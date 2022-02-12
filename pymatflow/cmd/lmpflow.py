#!/usr/bin/env python

import os
import sys
import argparse




def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="driver", title="subcommands", description="choose one and only one calculator")



    # --------------------------------------------------------------------------
    # LAMMPS
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("lmp", help="using lammps as calculator")

    subparser.add_argument("-r", "--runtype", type=int, default=0,
        choices=[0, 1, 2, 3, 4, 5, 6, 7],
        help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-cell; 3->hexagonal-cell; 4->tetragonal-cell; 5->neb; 6->vasp-phonon; 7->phonopy")

    subparser.add_argument("-d", "--directory", type=str, default="matflow-running",
        help="Directory for the running.")






    # --------------------------------------------------------------------------
    # LAMMPS data
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("data", help="lammps structure data operation")

    subparser.add_argument("-i", "--input", type=str,
        help="input structure file, canbe xxx.cif, xxx.xyz, xxx.xsd, xxx.cube")

    subparser.add_argument("-o", "--output", type=str, default="lammps.data",
        help="output lammps data file containing the structure information")
    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()

    # if no argument passed to matflow
    if len(sys.argv) == 1:
        # display help message when no args provided
        parser.print_help()
        sys.exit(1)


    if args.driver == "lmp":
        pass
    elif args.driver == "data":
        # will convert file type according to the suffix of the specified input
        # input structure
        if args.input.split(".")[-1] == "xyz":
            from pymatflow.structure.crystal import Crystal
            a = Crystal()
            a.from_xyz_file(args.input)
        elif args.input.split(".")[-1] == "cif":
            from pymatflow.structure.crystal import Crystal
            import pymatflow.third.aseio as aseio
            a = Crystal()
            a.cell, a.atoms = aseio.read_cif(args.input)
        elif args.input.split(".")[-1] == "xsd":
            from pymatflow.structure.crystal import Crystal
            import pymatflow.third.aseio as aseio
            a = Crystal()
            a.cell, a.atoms = aseio.read_xsd(args.input)
        elif args.input.split(".")[-1] == "xsf":
            from pymatflow.structure.crystal import Crystal
            import pymatflow.third.aseio as aseio
            a = Crystal()
            a.cell, a.atoms = aseio.read_xsf(args.input)
        elif os.path.basename(args.input) == "POSCAR" or os.path.basename(args.input) == "CONTCAR":
            from pymatflow.structure.crystal import Crystal
            import pymatflow.third.aseio as aseio
            a = Crystal()
            a.cell, a.atoms = aseio.read_poscar(args.input)

        # output structure
        import pymatflow.third.aseio
        aseio.write_lammps_data(cell=a.cell, atoms=a.atoms, filepath=args.output)

        print("=========================================================\n")
        print("              lmpflow data\n")
        print("---------------------------------------------------------\n")
        print("with the help from ase.io\n")



if __name__ == "__main__":
    main()
