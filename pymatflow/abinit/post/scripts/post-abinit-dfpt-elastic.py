#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import argparse
import copy
import datetime
import subprocess
import matplotlib.pyplot as plt

from pymatflow.abinit.post.dfpt import dfpt_elastic_anaddb_out



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously dfpt elastic running directory", type=str, default="tmp-abinit-dfpt-elastic")

    args = parser.parse_args()

    os.chdir(args.directory)
    elastic = dfpt_elastic_anaddb_out()
    elastic.get_info(file="anaddb.out")

    os.system("mkdir -p post-processing")
    os.chdir("post-processing")

    with open("elastic.md", 'w', encoding='utf-8') as fout:
        fout.write("# 弹性常数实验报告\n")
        fout.write("**Elastic Tensor (clamped ion)`%s`**\n" % elastic.run_info["elastic_tensor_clamped_ion_unit"])
        fout.write("```\n")
        for i in range(6):
            for j in range(6):
                fout.write("%.7f " % elastic.run_info["elastic_tensor_clamped_ion"][i][j])
            fout.write("\n")
        fout.write("```\n")
        fout.write("**Elastic Tensor (relaxed ion)`%s`**\n" % elastic.run_info["elastic_tensor_relaxed_ion_unit"])
        fout.write("```\n")
        for i in range(6):
            for j in range(6):
                fout.write("%.7f " % elastic.run_info["elastic_tensor_relaxed_ion"][i][j])
            fout.write("\n")
        fout.write("```\n")
        fout.write("**Compliance Tensor (clamped ion)`%s`**\n" % elastic.run_info["compliance_tensor_clamped_ion_unit"])
        fout.write("```\n")
        for i in range(6):
            for j in range(6):
                fout.write("%.7f " % elastic.run_info["compliance_tensor_clamped_ion"][i][j])
            fout.write("\n")
        fout.write("```\n")
        fout.write("**Compliance Tensor (relaxed ion)`%s`**\n" % elastic.run_info["compliance_tensor_relaxed_ion_unit"])
        fout.write("```\n")
        for i in range(6):
            for j in range(6):
                fout.write("%.7f " % elastic.run_info["compliance_tensor_relaxed_ion"][i][j])
            fout.write("\n")
        fout.write("```\n")
    os.chdir("../")
    os.chdir("../")


    # --------------------------------------------------------------------------
    # print information to the terminal
    # --------------------------------------------------------------------------
    print("=====================================================================\n")
    print("                           post-abinit-dfpt-elastic.py\n")
    print("---------------------------------------------------------------------\n")
    print("\n")
