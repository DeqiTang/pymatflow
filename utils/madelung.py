#!/usr/bin/env python

"""
calculate the madelung constant
"""
import sys
import argparse

from pymatflow.cmd.structflow import read_structure
from pymatflow.structure.crystal import crystal


"""
We should define a new xyz class that can read charge information from xyz file.
we set charge information(positive or negative) right after x y z coordinate of
each atom.
"""

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")

    args = parser.parse_args()

    # if no argument passed to matflow
    if len(sys.argv) == 1:
        # display help message when no args provided
        parser.print_help()
        sys.exit(1)
        
    a = read_structure(filepath=args.input)
    for i in range(1, 10):
        supercell = a.build_supercell([i, i, i])
        new_structure = crystal()
        new_structure.get_cell_atoms(cell=supercell["cell"], atoms=supercell["atoms"])
        # calculate the madelung constant
        # get center atom
        all_x = [atoms.x for atom in new_structure.atoms]
        all_y = [atoms.y for atom in new_structure.atoms]
        all_z = [atoms.z for atom in new_structure.atoms]
        
        

if __name__ == "__main__":
    main()