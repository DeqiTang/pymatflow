#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
import os
import shutil

from emuhelper.cp2k.base.glob import cp2k_glob
from emuhelper.cp2k.base.force_eval import cp2k_force_eval


"""
Usage:
    single_point_cp2k.py input.xyz
"""



if __name__ == "__main__":
    inp_name = "single-point.inp"
    glob = cp2k_glob()
    force_eval = cp2k_force_eval(sys.argv[1])

    if os.path.exists("tmp-single-point"):
        shutil.rmtree("tmp-single-point")
    os.mkdir("tmp-single-point")
    os.chdir("tmp-single-point")
    shutil.copyfile("../%s" % sys.argv[1], "./%s" % sys.argv[1])

    glob.params["RUN_TYPE"] = "ENERGY_FORCE"
    glob.to_global(inp_name)
    force_eval.to_force_eval(inp_name)
    os.system("cp2k.psmp -in %s | tee log.txt" % inp_name)
