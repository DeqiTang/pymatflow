#!/usr/bin/env python
# _*_ coding: utf-8 _*_


import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="intput xyz file", type=str)
parser.add_argument("-o", "--output", help="output xyz file", type=str)
parser.add_argument("--fix", help="list of fixed atoms, index start from 1", nargs='+', type=int)

args = parser.parse_args()

lines = None
with open(args.input, 'r') as fin:
    lines = fin.readlines()

for index in args.fix:
    lines[index+1] = lines[index+1].split("\n")[0] + " T T T\n"

with open(args.output, 'w') as fout:
    for line in lines:
        fout.write(line)
