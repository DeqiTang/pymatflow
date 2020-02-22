#!/usr/bin/env python
# _*_ coding: utf-8 _*_
import numpy as np

xyz_coord_file="lih-2x2x2.xyz"
siesta_coord_file="sieta_coord.txt"
species_total = {"H": 1, "Li": 2}

with open(siesta_coord_file, 'w') as fout:
  with open(xyz_coord_file, 'r') as f:
    i = 0
    for line in f:
      if i > 1:
        species_num = species_total[str(line.split()[0])]
        fout.write(str(line.split()[1]))
        fout.write("\t")
        fout.write(str(line.split()[2]))
        fout.write("\t")
        fout.write(str(line.split()[3]))
        fout.write("\t")
        fout.write(str(species_num))
        fout.write("\n")
      i += 1

