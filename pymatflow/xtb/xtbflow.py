#!/bin/env python

from pymatflow.xtb.md import md

def main():
    print("======================================================\n")
    print("               XTB command recommendation\n")
    print("------------------------------------------------------\n")
    print("Choose your type of calculation, supported are:\n")
    print("1) Molecular dynamics\n")
    type_of_calc = input()
    if type_of_calc == "1":
        md()


if __name__ == "__main__":
    main()
