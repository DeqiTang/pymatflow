#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import argparse
import copy
import datetime
import subprocess
import numpy as np
import matplotlib.pyplot as plt

from pymatflow.abinit.post.dfpt import dfpt_elastic_piezo_anaddb_out



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="previously dfpt elastic running directory", type=str, default="tmp-abinit-dfpt-elastic-piezo")

    args = parser.parse_args()

    os.chdir(args.directory)
    elastic = dfpt_elastic_piezo_anaddb_out()
    elastic.get_info(file="anaddb.out")

    os.system("mkdir -p post-processing")
    os.chdir("post-processing")

    with open("elastic-piezo.md", 'w', encoding='utf-8') as fout:
        fout.write("# 弹性张量-压电张量实验报告\n")
        fout.write("## 计算输出原始数据\n")
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

        fout.write("**Proper piezoelectric constants (clamped ion)`%s`**\n" % elastic.run_info["piezo_tensor_clamped_ion_unit"])
        fout.write("```\n")
        for i in range(6):
            for j in range(3):
                fout.write("%.8f " % elastic.run_info["piezo_tensor_clamped_ion"][i][j])
            fout.write("\n")
        fout.write("```\n")
        fout.write("**Proper piezoelectric constants (relaxed ion)`%s`**\n" % elastic.run_info["piezo_tensor_relaxed_ion_unit"])
        fout.write("```\n")
        for i in range(6):
            for j in range(3):
                fout.write("%.8f " % elastic.run_info["piezo_tensor_relaxed_ion"][i][j])
            fout.write("\n")
        fout.write("```\n")
        fout.write("\n")
        fout.write("**压电常数的计算如果ecut或者kmesh不够满足要求, 在anaddb输出压电张量时会附带警告:**\n")
        fout.write("```\n")
        fout.write("ddb_piezo : WARNING -\n")
        fout.write("Acoustic sum rule violation met : the eigenvalues of accoustic mode\n")
        fout.write("are too large at Gamma point\n")
        fout.write("Increase cutoff energy or k-points sampling.\n")
        fout.write("```\n")
        fout.write("你可以查看该文件\n")
        fout.write('\n')

        fout.write("## 数据处理\n")
        fout.write("如果这里anaddb分析出的Compliance tensor就是所谓的柔性张量那么我们就可以转换得到压电应变常数\n")
        fout.write("实际柔性张量就是弹性张量的逆矩阵, 我们可以检查一下上面的elastic tensor与compliance tensor矩阵乘积是否为单位矩阵\n")
        fout.write("测试clamped时elastic与compliance两个矩阵的乘积为:\n")
        multiply = np.mat(elastic.run_info["elastic_tensor_clamped_ion"]) * np.mat(elastic.run_info["compliance_tensor_clamped_ion"])
        fout.write("```\n")
        for i in range(6):
            for j in range(6):
                fout.write("%.7f " % multiply[i, j])
            fout.write("\n")
        fout.write("```\n")
        fout.write("的确是单位阵, 表明我们对其认识是正确的。\n")
        fout.write("**压电应变常数的计算**\n")
        fout.write("上面的压电常数时压电应力常数单位是 $C/m^(-2)$\n")
        fout.write("而压电应力常数与压电应变常数可以通过弹性柔顺系数进行转换:\n")
        fout.write("$[e] = [[d]^{t}[s]^{-1}]^{t}$\n")
        fout.write("$[d] = [[s][e]]^{t}$\n")
        fout.write("$[d] = [e]^{t}[s]^{t}$\n")
        fout.write("$$\n")

        e = np.mat(elastic.run_info["piezo_tensor_clamped_ion"])
        s = np.mat(elastic.run_info["compliance_tensor_clamped_ion"])
        d = np.transpose(s * e) * 1.0e4
        fout.write("计算所得压电应变常数为(clamped ion) $pC/N$:\n")
        fout.write("```\n")
        for i in range(3):
            for j in range(6):
                fout.write("%.8f " % d[i, j])
            fout.write("\n")
        fout.write("```\n")
        e = np.mat(elastic.run_info["piezo_tensor_relaxed_ion"])
        s = np.mat(elastic.run_info["compliance_tensor_relaxed_ion"])
        d = np.transpose(s * e) * 1.0e4
        fout.write("计算所得压电应变常数为(relaxed ion)$pC/N$:\n")
        fout.write("```\n")
        for i in range(3):
            for j in range(6):
                fout.write("%.8f " % d[i, j])
            fout.write("\n")
        fout.write("```\n")


    os.chdir("../")
    os.chdir("../")


    # --------------------------------------------------------------------------
    # print information to the terminal
    # --------------------------------------------------------------------------
    print("=====================================================================\n")
    print("                           post-abinit-dfpt-elastic-piezo.py\n")
    print("---------------------------------------------------------------------\n")
    print("\n")
