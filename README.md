# emuhelper

With scripts to provide preparing and/or post-processing for your Emulations(DFT, MD, FEM).


## TODO
### CP2K
* implementing Linear Scaling SCF [OK]
* implementing VDW correction
* implementing setting of kpoints for band structure and dos calculation
* implementing neb calculation script module
* implementing `MULTI_FORCE_EVALS
* implementing classic mocular dynamics
* implementing QMMM
* implementing Metadynamics
* implementing VIBRATIONAL ANALYSIS
* implementing Excited States calculation using TDDFPT
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

## About
其实有GUI的构建工具比如Avogadro、Pymol等可以用于构建一些模型。但是我想要构建一个自动化程度比较高的工具来实现简化模拟的准备工作。
