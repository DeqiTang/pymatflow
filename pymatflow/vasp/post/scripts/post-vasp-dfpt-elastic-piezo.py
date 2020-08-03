#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import argparse
import copy
import datetime
import subprocess
import numpy as np
import matplotlib.pyplot as plt

"""
Piezoelastic Strain Tensor: IBIRION = 8, LEPSILON = T
Elastic tensor: IBIRION = 6, LEPSILON = T, NFREE = 4, ISIF = 3, to get the TOTAL ELASTIC MODULI, ISIF = 3 and NFREE = 4 is needed.
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outcars", help="OUTCAR for piezoelectric stress tensor and elastic tensor calculation respectively", type=str, nargs=2, required=True)

    args = parser.parse_args()


    # extract e_from_electrons, e_from_ions and c_elastic: converted to numpy array
    
    # OUTCAR to extract piezoelectric stress tensor
    with open(os.path.join(args.outcars[0]), "r") as fin:
        outcar_lines = fin.readlines()
        
    e_from_electrons = []
    e_from_ions = []
    for i in range(len(outcar_lines)):
        if len(outcar_lines[i].split()) == 0:
            continue
        if outcar_lines[i].split()[0] == "PIEZOELECTRIC" and outcar_lines[i].split()[1] == "TENSOR" and outcar_lines[i].split()[2] == "for" and outcar_lines[i].split()[8].split("\n")[0] == "(C/m^2)":
            for j in range(i+3, i+6, 1):
                tmp = []
                for k in range(6):
                    tmp.append(float(outcar_lines[j].split()[k+1]))
                e_from_electrons.append(tmp)
        if outcar_lines[i].split()[0] == "PIEZOELECTRIC" and outcar_lines[i].split()[1] == "TENSOR" and outcar_lines[i].split()[2] == "IONIC" and outcar_lines[i].split()[10].split("\n")[0] == "(C/m^2)":
            for j in range(i+3, i+6, 1):
                tmp = []
                for k in range(6):
                    tmp.append(float(outcar_lines[j].split()[k+1]))
                e_from_ions.append(tmp)                 
    e_from_electrons = np.array(e_from_electrons)              
    e_from_ions = np.array(e_from_ions)              

    # OUTCAR to extract elastic tensor
    with open(os.path.join(args.outcars[1]), "r") as fin:
        outcar_lines = fin.readlines()
        
    c_elastic = []
    for i in range(len(outcar_lines)):
        if len(outcar_lines[i].split()) == 0:
            continue
        if outcar_lines[i].split()[0] == "TOTAL" and outcar_lines[i].split()[1] == "ELASTIC" and outcar_lines[i].split()[2] == "MODULI" and outcar_lines[i].split()[3].split("\n")[0] == "(kBar)":
            for j in range(i+3, i+9, 1):
                tmp = []
                for k in range(6):
                    tmp.append(float(outcar_lines[j].split()[k+1]))
                c_elastic.append(tmp)
    c_elastic = np.array(c_elastic)      

    # e_from_electrons:
    # e_from_ions:
    # e_total = e_from_electrons + e_from_ions: C/m^2
    # c_elastic: kBar = 10^8 N/m^2
    # d_ = e_total * c_elastic^-1
    e_total = e_from_electrons + e_from_ions
    d_piezo_strain = np.dot(e_total, np.linalg.inv(c_elastic * 1.0e+8)) # C/N
    d_piezo_strain_pc_n = np.dot(e_total, np.linalg.inv(c_elastic * 1.0e+8)) * 1.0e+12 # pC/N
        
    
    print("*********************************************************\n")
    print("       Calculated piezoelectric strain constant\n")
    print("---------------------------------------------------------\n")
    print("Piezoelectric strain tensor (C/N)\n")
    print("%f %f %f %f %f %f\n" % (d_piezo_strain[0, 0], d_piezo_strain[0, 1], d_piezo_strain[0, 2], d_piezo_strain[0, 3], d_piezo_strain[0, 4], d_piezo_strain[0, 5]))
    print("%f %f %f %f %f %f\n" % (d_piezo_strain[1, 0], d_piezo_strain[1, 1], d_piezo_strain[1, 2], d_piezo_strain[1, 3], d_piezo_strain[1, 4], d_piezo_strain[1, 5]))
    print("%f %f %f %f %f %f\n" % (d_piezo_strain[2, 0], d_piezo_strain[2, 1], d_piezo_strain[2, 2], d_piezo_strain[2, 3], d_piezo_strain[2, 4], d_piezo_strain[2, 5]))
    print("Piezoelectric strain tensor (pC/N)\n")
    print("%f %f %f %f %f %f\n" % (d_piezo_strain_pc_n[0, 0], d_piezo_strain_pc_n[0, 1], d_piezo_strain_pc_n[0, 2], d_piezo_strain_pc_n[0, 3], d_piezo_strain_pc_n[0, 4], d_piezo_strain_pc_n[0, 5]))
    print("%f %f %f %f %f %f\n" % (d_piezo_strain_pc_n[1, 0], d_piezo_strain_pc_n[1, 1], d_piezo_strain_pc_n[1, 2], d_piezo_strain_pc_n[1, 3], d_piezo_strain_pc_n[1, 4], d_piezo_strain_pc_n[1, 5]))
    print("%f %f %f %f %f %f\n" % (d_piezo_strain_pc_n[2, 0], d_piezo_strain_pc_n[2, 1], d_piezo_strain_pc_n[2, 2], d_piezo_strain_pc_n[2, 3], d_piezo_strain_pc_n[2, 4], d_piezo_strain_pc_n[2, 5]))
    
