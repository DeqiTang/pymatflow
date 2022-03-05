#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import re
import subprocess
from glob import glob

import setuptools
from setuptools import setup
from setuptools import find_packages, Extension
from setuptools import Command
from setuptools.command.build_ext import build_ext

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')



setup(
    name = "pymatflow",
    version = '0.1.1a4',
    ## python3 setup.py build sdist bdist_wheel
    ## twine upload dist/*
    keywords = ("Ab intio ,DFT, workflow, input generation"),
    description = "An emulation assistant, input generation and manage for DFT programs",
    license = "MIT",
    url = "https://gitlab.com/deqitang/pymatflow",
    author_email = "pymatflow@163.com",
    packages = find_packages(),
    #data_files
    include_package_data = True,
    platforms = "any",
    python_requires = ">=3.0",
    install_requires = ["atomsciflow", "ase", "numpy", "scipy", "matplotlib"],
    cmdclass = {
        'clean': CleanCommand,
        'build_ext': build_ext
    },    
    ext_modules=[
    ],
    # ------------------------------------------
    scripts = [
        "pymatflow/scripts/xyz-build-supercell.py",
        "pymatflow/scripts/cluster_sphere.py",
        "pymatflow/scripts/xyz-fix-atoms.py",
        "pymatflow/scripts/seekpath-flow.py",
        "pymatflow/scripts/cif-to-pdb.py",
        "pymatflow/scripts/pdb-to-cif.py",
        "pymatflow/scripts/xyz-modified-to-crystal.py",
        "pymatflow/scripts/kpath-seekpath.py",
        "pymatflow/scripts/kpath-cms-2010.py",
        "pymatflow/scripts/kpoints-mp-suggest.py",
        "pymatflow/scripts/pot-from-xyz-modified.py",
        "pymatflow/scripts/vasp-potcar-from-xyz.py",
        "pymatflow/scripts/xyz-modified-to-poscar.py",
        "pymatflow/scripts/llhpc-env.py",
        "pymatflow/scripts/nebmake.py",
        "pymatflow/scripts/contour-xyz.py",
        "pymatflow/scripts/contour-xyz-diff.py",
        "pymatflow/scripts/ldos-xyz.py",  
        "pymatflow/scripts/xsd-color-bond.py",
        "pymatflow/scripts/stm-vasp.py",
        "pymatflow/scripts/chg-vasp.py",
        "pymatflow/scripts/cube-handle.py",
        "pymatflow/scripts/chg-diff-1d.py",
        "pymatflow/scripts/chgflow.py",
        "pymatflow/scripts/cube-diff-1d.py",
        "pymatflow/scripts/cube-1d.py",
        "pymatflow/scripts/veusz-flow.py",
        "pymatflow/remote/scripts/thq.py",
        "pymatflow/remote/scripts/thcancel.py",
        "pymatflow/remote/scripts/thpull.py",
        "pymatflow/remote/scripts/thcmd.py",
        "pymatflow/remote/scripts/sz-pull.py",
        "pymatflow/remote/scripts/sz-cmd.py",
        "pymatflow/remote/scripts/sz-del.py",
        "pymatflow/remote/scripts/sz-q.py",
        "pymatflow/remote/scripts/threport.py",
        "pymatflow/flow/scripts/flow-run.py",
        "pymatflow/flow/scripts/flow-pes-static-qe.py",
        "pymatflow/flow/scripts/flow-pes-relax-qe.py",
        "pymatflow/flow/scripts/flow-pes-relax-cp2k.py",
        "pymatflow/flow/scripts/pes-to-img.py",
        "pymatflow/vasp/scripts/vasp-time.py",
        "pymatflow/vasp/scripts/vasp-grep-energy.py",
        "pymatflow/vasp/scripts/xdatcar-to-xyz.py",
        "pymatflow/vasp/post/scripts/vasp-tensor-piezo.py",
        "pymatflow/octopus/scripts/octflow.py",
        ],
    entry_points = {
        'console_scripts': [
            'matflow = pymatflow.cmd.matflow:main',
            'mflow = pymatflow.cmd.matflow:main',
            'postflow = pymatflow.cmd.postflow:main',
            'pflow = pymatflow.cmd.postflow:main',
            'structflow = pymatflow.cmd.structflow:main',
            'sflow = pymatflow.cmd.structflow:main',
            'lmpflow = pymatflow.cmd.lmpflow:main',
            'lflow = pymatflow.cmd.lmpflow:main',
            "psoflow = pymatflow.flow.calypso.psoflow:main",
            "octoflow = pymatflow.octopus.octoflow:main",
            "excitingflow = pymatflow.exciting.excitingflow:main",
            "elkflow = pymatflow.elk.elkflow:main",
            "xtbflow = pymatflow.xtb.xtbflow:main",
        ]
    },
)

