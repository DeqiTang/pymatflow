# emuhelper

With scripts to provide preparing and/or post-processing for your Emulations(DFT, MD, FEM).

## Goals
* powerful and stable calculator supports for popular Ab initio programs
* fully featured results analysis
* user friendly structure builder
* automated workflow for high throughput calculation[including database handling]
  * automatical choice of appropriate parameters
  * automatical correction of error during running
* statistical analysis for properties, basic machine learning module

## TODO
### CP2K
* implementing DFTB calculation
* implementing Quantum Espresso type calculation by using SIRIUS in cp2k
* implementing Linear Scaling SCF [OK]
* implementing VDW correction [OK]
* implementing setting of kpoints for band structure and dos calculation
* implementing neb calculation script module [OK]
* implementing `MULTI_FORCE_EVALS
    * we can run optimization, neb, md based on MULTI_FORCE_EVALS
    * Reference: [Nanostructures and adsorption on metallic surfaces](https://www.cp2k.org/exercises:2015_cecam_tutorial:neb)
* implementing classic mocular dynamics
* implementing QMMM
* implementing Metadynamics
* implementing Monte Carlo simulation
* implementing EHRENFEST_DYN
* implementing LINEAR_RESPONSE [OK]
* implementing NEGF
* implementing MC
* implementing PINT
* implementing SPECTRA
* implementing MP2 method
* implementing VIBRATIONAL ANALYSIS [OK]
* implementing Excited States calculation using TDDFPT
* implementing []TDDFT for XAS](https://www.cp2k.org/exercises:2019_conexs_newcastle:index)
* implementing Transport calculations with NEGF
* implementing RESP Charges analysis [OK]
* develop post processing module
### Quantum Espresso
* improving Phonon calculation using ph.x
* implementing kpoints setting for band structure calculation
* develop post processing module
### SIESTA
* improving code base structure
* develop post processing module
### Abinit
* considering reconstructing the code for multi-dataset using
* develop post processing module
### ORCA
* build the base structure
* develop post processing module
### dalton
* build the base structure
* develop post processing module



## Features
* Automated surface cut: not implemented yet
* Slab Builder: not implemented yet
* Nanotube Builder: not implemented yet
* Cluster Builder: not implemented yet
* Grain Boundary Builder: not implemented yet
* Polycrystalline Builder(full atom): not implemented yet
* Automated Hybridization: not implemented yet
* Domain Boundary Design: not implemented yet
* Super Lattice: not implemented yet
* Transformation Matrix: not implemented yet
* Support over 1000,000 atoms: not sure now
* Auto determine the cell parameters for the xyz structure: not implemented yet
* More simple scripts to working with remote server

## About the th\*.py scripts
in order to use th\*.py scripts and the function of --auto in running scripts(like qe-relax.py), you have to prepare a server.conf file in ~/.emuhelper/server.conf, in format like this:
```
[server]
ip = x.x.x.x
user = xxx
serverdir = /xx/xx/xx # like /THFS/home/xxx/
password = xxx
```

## About
其实有GUI的构建工具比如Avogadro、Pymol等可以用于构建一些模型。但是我想要构建一个自动化程度比较高的工具来实现简化模拟的准备工作。
