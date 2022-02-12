#!/usr/bin/env python
import numpy as np
import argparse
from pymatflow.cmd.structflow import read_structure


def suggest(x, limit=17):
    """
    :params x: -> float
        the length of a or b or c of the crystal.

    :params limit: -> float
        the minimum limit for product of x and k.
        So k is the minimum positive integer that 
        satisfy: x * k >= limit
    Note:
        Rule for k grid suggestion:
        for x = a|b|c, k is the minimum positive integer that
        satisfy: x * k >= limit.
    """    
    k = 1
    while k * x < limit:
        k += 1
    return k

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, required=True,
            help="input structure file")
    
    parser.add_argument("-l", "--limit", type=float, default=17,
            help="k is the minimum that obey the rule: x * k >= limit, it is suggested that limit should be set from 15 to 30, default is 17")
    # ===============================================================================
    args = parser.parse_args()
    
    structure = read_structure(filepath=args.input)
    a = np.linalg.norm(structure.cell[0])
    b = np.linalg.norm(structure.cell[1])
    c = np.linalg.norm(structure.cell[2])

    print("=====================================================\n")
    print("        Monkhorst-Pack kpoints suggestion\n")
    print("-----------------------------------------------------\n")
    print("Rule for k grid suggestion:\n")
    print("for x = a|b|c, k is the minimum positive integer that\n")
    print("satisfy: x * k >= limit. and your value for limit is %f\n" % args.limit)
    print("The suggested kpoints grid is:\n")
    print("%d %d %d" % (
        suggest(a, limit=args.limit), 
        suggest(b, limit=args.limit), 
        suggest(c, limit=args.limit)
        ))

if __name__ == "__main__":
    main()
