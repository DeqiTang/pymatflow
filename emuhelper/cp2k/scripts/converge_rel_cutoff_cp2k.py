#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from emuhelper.cp2k.base.glob import cp2k_glob
from emuhelper.cp2k.base.force_eval import cp2k_force_eval

"""
Usage:
    python converge_rel_cutoff.py xxx.xyz rel_cutoff_min rel_cutoff_max rel_cutoff_step cutoff
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.
"""


rel_cutoff_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
rel_cutoff_max = int(sys.argv[3])
rel_cutoff_step = int(sys.argv[4])
cutoff = int(sys.argv[5])


if __name__ == '__main__':
    glob = cp2k_glob()
    force_eval = cp2k_force_eval(sys.argv[1])
    glob.params["RUN_TYPE"] = "ENERGY_FORCE"

    if os.path.exists("./tmp-rel-cutoff"):
        shutil.rmtree("./tmp-rel-cutoff")
    os.mkdir("./tmp-rel-cutoff")
    os.chdir("./tmp-rel-cutoff")
    shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])

    n_test = int((rel_cutoff_max - rel_cutoff_min) / rel_cutoff_step)
    for i in range(n_test + 1):
        rel_cutoff = int(rel_cutoff_min + i * rel_cutoff_step)
        inp_name = "test-rel-cutoff-%d.inp" % rel_cutoff
        glob.to_global(inp_name)
        force_eval.dft.mgrid.params["CUTOFF"] = cutoff
        force_eval.dft.mgrid.params["REL_CUTOFF"] = rel_cutoff
        force_eval.to_force_eval(inp_name)
    # run the simulation
    for i in range(n_test + 1):
        rel_cutoff = int(rel_cutoff_min + i * rel_cutoff_step)
        inp_name = "test-rel-cutoff-%d.inp" % rel_cutoff
        out_f_name = "test-rel-cutoff-%d.out" % rel_cutoff
        os.system("cp2k.psmp -in %s | tee %s" % (inp_name, out_f_name))
    # analyse the result
    for i in range(n_test + 1):
        rel_cutoff = int(rel_cutoff_min + i * rel_cutoff_step)
        out_f_name = "test-rel-cutoff-%d.out" % rel_cutoff
        os.system("cat %s | grep 'Total energy:' >> energy-rel-cutoff.data" % out_f_name)

    rel_cut = [rel_cutoff_min + i * rel_cutoff_step for i in range(n_test + 1)]
    energy = []
    with open("energy-rel-cutoff.data", 'r') as fin:
        for line in fin:
            energy.append(float(line.split()[2]))

    import matplotlib.pyplot as plt

    plt.plot(rel_cut, energy)
    plt.show()
