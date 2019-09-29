#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import re

"""
usage:
    post-single-point-cp2k.py xxx.out
    xxx.out is the output file of cp2k single point running
"""

def get_final_energy(fin):
    # fin: a file stream for reading (the output of static scf running)
    #energy_string = " ENERGY| Total FORCE_EVAL ( QS  ) energy (a.u.):*"
    energy_string = "\sENERGY| Total FORCE_EVAL ( QS  ) energy (a.u.):*"
    for line in fin:
        match = re.match(energy_string, line)
        if match is not None:
            return float(match.string.split()[8])

if __name__ == '__main__':
    with open(sys.argv[1], 'r') as fin:
        final_energy = get_final_energy(fin)
    print("The final free energy extrapolated for TS -> 0 is: %s Hartree = %s eV" % (str(final_energy), str(27.2114 * final_energy)))
