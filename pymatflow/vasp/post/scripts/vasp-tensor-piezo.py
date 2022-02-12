#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
import csv
import copy
import datetime
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pymatflow.cmd.structflow import read_structure

"""
Piezoelastic Strain Tensor: IBIRION = 8, LEPSILON = T
Elastic tensor: IBIRION = 6, LEPSILON = T, NFREE = 4, ISIF = 3, to get the TOTAL ELASTIC MODULI, ISIF = 3 and NFREE = 4 is needed.
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outcars", help="OUTCAR for piezoelectric stress tensor and elastic tensor calculation respectively", type=str, nargs=2, required=True)

    parser.add_argument("--poscar", type=str, required=True, help="POSCAR of the structure")
    
    parser.add_argument("--output-csv", type=str, default="./piezo_elastic_data.csv",
            help="specify the path for the csv file to store the results")

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
        #if outcar_lines[i].split()[0] == "PIEZOELECTRIC" and outcar_lines[i].split()[1] == "TENSOR" and outcar_lines[i].split()[2] == "for" and outcar_lines[i].split()[8].split("\n")[0] == "(C/m^2)":
        possible1 = "PIEZOELECTRIC TENSOR  for field in x, y, z        (C/m^2)"
        possible2 = "PIEZOELECTRIC TENSOR (including local field effects)  for field in x, y, z        (C/m^2)"
        # different version of VASP or different method may generate different output
        # so possible1 and possible2 should be dealt with at the same time.
        if possible1 in outcar_lines[i] or possible2 in outcar_lines[i]:
            for j in range(i+3, i+6, 1):
                tmp = []
                for k in range(6):
                    tmp.append(float(outcar_lines[j].split()[k+1]))
                e_from_electrons.append(tmp)
        #if outcar_lines[i].split()[0] == "PIEZOELECTRIC" and outcar_lines[i].split()[1] == "TENSOR" and outcar_lines[i].split()[2] == "IONIC" and outcar_lines[i].split()[10].split("\n")[0] == "(C/m^2)":
        if "PIEZOELECTRIC TENSOR IONIC CONTR  for field in x, y, z        (C/m^2)" in outcar_lines[i]:
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
        #if outcar_lines[i].split()[0] == "TOTAL" and outcar_lines[i].split()[1] == "ELASTIC" and outcar_lines[i].split()[2] == "MODULI" and outcar_lines[i].split()[3].split("\n")[0] == "(kBar)":
        if "TOTAL ELASTIC MODULI (kBar)" in outcar_lines[i]:
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
    e_total_pc_m2 = e_total * 1.0e+12 # pC/m^2
    
    print("*********************************************************\n")
    print("       Calculated piezoelectric strain constant\n")
    print("---------------------------------------------------------\n")
    print("Piezoelectric stress tensor (C/m^2)\n")
    print("%f %f %f %f %f %f\n" % (e_total[0, 0], e_total[0, 1], e_total[0, 2], e_total[0, 3], e_total[0, 4], e_total[0, 5]))
    print("%f %f %f %f %f %f\n" % (e_total[1, 0], e_total[1, 1], e_total[1, 2], e_total[1, 3], e_total[1, 4], e_total[1, 5]))
    print("%f %f %f %f %f %f\n" % (e_total[2, 0], e_total[2, 1], e_total[2, 2], e_total[2, 3], e_total[2, 4], e_total[2, 5]))
    
    print("Piezoelectric stress tensor (pC/m^2)\n")
    print("%f %f %f %f %f %f\n" % (e_total_pc_m2[0, 0], e_total_pc_m2[0, 1], e_total_pc_m2[0, 2], e_total_pc_m2[0, 3], e_total_pc_m2[0, 4], e_total_pc_m2[0, 5]))
    print("%f %f %f %f %f %f\n" % (e_total_pc_m2[1, 0], e_total_pc_m2[1, 1], e_total_pc_m2[1, 2], e_total_pc_m2[1, 3], e_total_pc_m2[1, 4], e_total_pc_m2[1, 5]))
    print("%f %f %f %f %f %f\n" % (e_total_pc_m2[2, 0], e_total_pc_m2[2, 1], e_total_pc_m2[2, 2], e_total_pc_m2[2, 3], e_total_pc_m2[2, 4], e_total_pc_m2[2, 5]))
    
    print("Piezoelectric strain tensor (C/N)\n")
    print("%f %f %f %f %f %f\n" % (d_piezo_strain[0, 0], d_piezo_strain[0, 1], d_piezo_strain[0, 2], d_piezo_strain[0, 3], d_piezo_strain[0, 4], d_piezo_strain[0, 5]))
    print("%f %f %f %f %f %f\n" % (d_piezo_strain[1, 0], d_piezo_strain[1, 1], d_piezo_strain[1, 2], d_piezo_strain[1, 3], d_piezo_strain[1, 4], d_piezo_strain[1, 5]))
    print("%f %f %f %f %f %f\n" % (d_piezo_strain[2, 0], d_piezo_strain[2, 1], d_piezo_strain[2, 2], d_piezo_strain[2, 3], d_piezo_strain[2, 4], d_piezo_strain[2, 5]))
    print("Piezoelectric strain tensor (pC/N)\n")
    print("%f %f %f %f %f %f\n" % (d_piezo_strain_pc_n[0, 0], d_piezo_strain_pc_n[0, 1], d_piezo_strain_pc_n[0, 2], d_piezo_strain_pc_n[0, 3], d_piezo_strain_pc_n[0, 4], d_piezo_strain_pc_n[0, 5]))
    print("%f %f %f %f %f %f\n" % (d_piezo_strain_pc_n[1, 0], d_piezo_strain_pc_n[1, 1], d_piezo_strain_pc_n[1, 2], d_piezo_strain_pc_n[1, 3], d_piezo_strain_pc_n[1, 4], d_piezo_strain_pc_n[1, 5]))
    print("%f %f %f %f %f %f\n" % (d_piezo_strain_pc_n[2, 0], d_piezo_strain_pc_n[2, 1], d_piezo_strain_pc_n[2, 2], d_piezo_strain_pc_n[2, 3], d_piezo_strain_pc_n[2, 4], d_piezo_strain_pc_n[2, 5]))
    
    print("Total Elastic tensor (kBar)\n")
    print("%f %f %f %f %f %f\n" % (c_elastic[0, 0], c_elastic[0, 1], c_elastic[0, 2], c_elastic[0, 3], c_elastic[0, 4], c_elastic[0, 5]))
    print("%f %f %f %f %f %f\n" % (c_elastic[1, 0], c_elastic[1, 1], c_elastic[1, 2], c_elastic[1, 3], c_elastic[1, 4], c_elastic[1, 5]))
    print("%f %f %f %f %f %f\n" % (c_elastic[2, 0], c_elastic[2, 1], c_elastic[2, 2], c_elastic[2, 3], c_elastic[2, 4], c_elastic[2, 5]))
    print("%f %f %f %f %f %f\n" % (c_elastic[3, 0], c_elastic[3, 1], c_elastic[3, 2], c_elastic[3, 3], c_elastic[3, 4], c_elastic[3, 5]))
    print("%f %f %f %f %f %f\n" % (c_elastic[4, 0], c_elastic[4, 1], c_elastic[4, 2], c_elastic[4, 3], c_elastic[4, 4], c_elastic[4, 5]))
    print("%f %f %f %f %f %f\n" % (c_elastic[5, 0], c_elastic[5, 1], c_elastic[5, 2], c_elastic[5, 3], c_elastic[5, 4], c_elastic[5, 5]))
    print("Total Elastic tensor (N/m^2)\n")
    c_elastic_n_m2 = c_elastic * 1.0e+8
    print("%f %f %f %f %f %f\n" % (c_elastic_n_m2[0, 0], c_elastic_n_m2[0, 1], c_elastic_n_m2[0, 2], c_elastic_n_m2[0, 3], c_elastic_n_m2[0, 4], c_elastic_n_m2[0, 5]))
    print("%f %f %f %f %f %f\n" % (c_elastic_n_m2[1, 0], c_elastic_n_m2[1, 1], c_elastic_n_m2[1, 2], c_elastic_n_m2[1, 3], c_elastic_n_m2[1, 4], c_elastic_n_m2[1, 5]))
    print("%f %f %f %f %f %f\n" % (c_elastic_n_m2[2, 0], c_elastic_n_m2[2, 1], c_elastic_n_m2[2, 2], c_elastic_n_m2[2, 3], c_elastic_n_m2[2, 4], c_elastic_n_m2[2, 5]))
    print("%f %f %f %f %f %f\n" % (c_elastic_n_m2[3, 0], c_elastic_n_m2[3, 1], c_elastic_n_m2[3, 2], c_elastic_n_m2[3, 3], c_elastic_n_m2[3, 4], c_elastic_n_m2[3, 5]))
    print("%f %f %f %f %f %f\n" % (c_elastic_n_m2[4, 0], c_elastic_n_m2[4, 1], c_elastic_n_m2[4, 2], c_elastic_n_m2[4, 3], c_elastic_n_m2[4, 4], c_elastic_n_m2[4, 5]))
    print("%f %f %f %f %f %f\n" % (c_elastic_n_m2[5, 0], c_elastic_n_m2[5, 1], c_elastic_n_m2[5, 2], c_elastic_n_m2[5, 3], c_elastic_n_m2[5, 4], c_elastic_n_m2[5, 5]))
    
    print("Piezoelectric stress tensor of 2D materials (z as surface direction)\n")
    print("Piezoelectric stess tensor in 2D (pC/m)\n")
    print("e_ij(2D)\n")
    print("e11 e12 e16\n")
    print("e21 e22 e26\n")
    print("e31 e32 e36\n")    
    structure = read_structure(filepath=args.poscar)
    c = np.linalg.norm(np.array(structure.cell[2]))
    e_total_2d_pc_m = []
    for i in [1, 2, 3]:
        row = []
        for j in [1, 2, 6]:
            row.append(e_total[i-1, j-1])
        e_total_2d_pc_m.append(row)
    e_total_2d_pc_m = np.array(e_total_2d_pc_m) * 1.0e+12 * c * 1.0e-10 # pC/m
    print("%f %f %f\n" % (e_total_2d_pc_m[0, 0], e_total_2d_pc_m[0, 1], e_total_2d_pc_m[0, 2]))
    print("%f %f %f\n" % (e_total_2d_pc_m[1, 0], e_total_2d_pc_m[1, 1], e_total_2d_pc_m[1, 2]))
    print("%f %f %f\n" % (e_total_2d_pc_m[2, 0], e_total_2d_pc_m[2, 1], e_total_2d_pc_m[2, 2]))
    print("\n")
    print("e_ij(2D): absolute value\n")
    print("| |e11| |e12| |e16| |\n")
    print("| |e21| |e22| |e26| |\n")
    print("| |e31| |e32| |e36| |\n")
    e_total_2d_pc_m_abs = np.abs(e_total_2d_pc_m) # we use absolute value of eij in 2d condition
    print("%f %f %f\n" % (e_total_2d_pc_m_abs[0, 0], e_total_2d_pc_m_abs[0, 1], e_total_2d_pc_m_abs[0, 2]))
    print("%f %f %f\n" % (e_total_2d_pc_m_abs[1, 0], e_total_2d_pc_m_abs[1, 1], e_total_2d_pc_m_abs[1, 2]))
    print("%f %f %f\n" % (e_total_2d_pc_m_abs[2, 0], e_total_2d_pc_m_abs[2, 1], e_total_2d_pc_m_abs[2, 2]))
    
    print("Elastic tensor in 2D (N/m)\n")
    print("C_ij(2D)\n")
    print("C11 C12 C16\n")
    print("C21 C22 C26\n")
    print("C61 C62 C66\n")
    # c_elastic: kBar = 1.0e+8 N/m^2
    c_elastic_2d_n_m = []
    for i in [1, 2, 6]:
        row = []
        for j in [1, 2, 6]:
            row.append(c_elastic[i-1, j-1])
        c_elastic_2d_n_m.append(row)
    c_elastic_2d_n_m = np.array(c_elastic_2d_n_m) *1.0e+8 * c * 1.0e-10 # N/m
    print("%f %f %f\n" % (c_elastic_2d_n_m[0, 0], c_elastic_2d_n_m[0, 1], c_elastic_2d_n_m[0, 2]))
    print("%f %f %f\n" % (c_elastic_2d_n_m[1, 0], c_elastic_2d_n_m[1, 1], c_elastic_2d_n_m[1, 2]))
    print("%f %f %f\n" % (c_elastic_2d_n_m[2, 0], c_elastic_2d_n_m[2, 1], c_elastic_2d_n_m[2, 2]))
    
    print("Piezoelectric strain tensor in 2D (pC/N)\n")
    #d_piezo_strain_2d_pc_n = np.dot(e_total_2d_pc_m, np.linalg.inv(c_elastic_2d_n_m)) # pC/N
    #print("%f %f %f\n" % (d_piezo_strain_2d_pc_n[0, 0], d_piezo_strain_2d_pc_n[0, 1], d_piezo_strain_2d_pc_n[0, 2]))
    #print("%f %f %f\n" % (d_piezo_strain_2d_pc_n[1, 0], d_piezo_strain_2d_pc_n[1, 1], d_piezo_strain_2d_pc_n[1, 2]))
    #print("%f %f %f\n" % (d_piezo_strain_2d_pc_n[2, 0], d_piezo_strain_2d_pc_n[2, 1], d_piezo_strain_2d_pc_n[2, 2]))
    d_piezo_strain_2d_pc_n = np.dot(e_total_2d_pc_m_abs, np.linalg.inv(c_elastic_2d_n_m)) # pC/N
    print("%f %f %f\n" % (d_piezo_strain_2d_pc_n[0, 0], d_piezo_strain_2d_pc_n[0, 1], d_piezo_strain_2d_pc_n[0, 2]))
    print("%f %f %f\n" % (d_piezo_strain_2d_pc_n[1, 0], d_piezo_strain_2d_pc_n[1, 1], d_piezo_strain_2d_pc_n[1, 2]))
    print("%f %f %f\n" % (d_piezo_strain_2d_pc_n[2, 0], d_piezo_strain_2d_pc_n[2, 1], d_piezo_strain_2d_pc_n[2, 2]))
    
    
    # output data to csv file
    with open(args.output_csv, "w") as fout:
        csv_writer = csv.writer(fout)
        csv_writer.writerow(["Piezoelectric stress tensor (C/m^2)"])
        csv_writer.writerows(e_total.tolist())
        
        csv_writer.writerow(["Piezoelectric stress tensor (pC/m^2)"])
        csv_writer.writerows(e_total_pc_m2.tolist())
        
        csv_writer.writerow(["Piezoelectric strain tensor (C/N)"])
        csv_writer.writerows(d_piezo_strain.tolist())

        csv_writer.writerow(["Piezoelectric strain tensor (pC/N)"])
        csv_writer.writerows(d_piezo_strain_pc_n.tolist())

        csv_writer.writerow(["Total Elastic tensor (kBar)"])
        csv_writer.writerows(c_elastic.tolist())

        csv_writer.writerow(["Total Elastic tensor (N/m^2)"])
        csv_writer.writerows(c_elastic_n_m2.tolist())
    
        csv_writer.writerow(["Piezoelectric stress tensor of 2D materials (z as surface direction)"])
        csv_writer.writerow(["Piezoelectric stess tensor in 2D (pC/m)"])
        csv_writer.writerow(["e_ij(2D)"])
        csv_writer.writerow(["e11 e12 e16"])
        csv_writer.writerow(["e21 e22 e26"])
        csv_writer.writerow(["e31 e32 e36"])
        csv_writer.writerows(e_total_2d_pc_m.tolist())    
    
        csv_writer.writerow(["e_ij(2D): absolute value"])
        csv_writer.writerow(["| |e11| |e12| |e16| |"])
        csv_writer.writerow(["| |e21| |e22| |e26| |"])
        csv_writer.writerow(["| |e31| |e32| |e36| |"])
        csv_writer.writerows(e_total_2d_pc_m_abs.tolist())
    
    
        csv_writer.writerow(["Elastic tensor in 2D (N/m)"])
        csv_writer.writerow(["C_ij(2D)"])
        csv_writer.writerow(["C11 C12 C16"])
        csv_writer.writerow(["C21 C22 C26"])
        csv_writer.writerow(["C61 C62 C66"])
        csv_writer.writerows(c_elastic_2d_n_m.tolist())
        
        csv_writer.writerow(["Piezoelectric strain tensor in 2D (pC/N)"])
        csv_writer.writerows(d_piezo_strain_2d_pc_n.tolist())