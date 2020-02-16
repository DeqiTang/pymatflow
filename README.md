# Pymatflow

Pymatflow is a workflow manager and highthroughput tool for DFT modeling of materials. It can automatically preparing input files for popular open source DFT programs like Quantum ESPRESSO, CP2K, Abinit, and SIESTA. And the development for the post processing is on the way.

## Install
To install the current stable release, just use pip to install Pymatflow from PyPI:
```
$ pip install pymatflow
```
to install the current develpment version, you can clone the project and manually install it after chdir to the project root directory:
```
$ git clone https://gitlab.com/deqitang/pymatflow.git
$ cd pymatflow
$ python setup.py install
```
or directly:
```
$ pip install git+https://gitlab.com/deqitang/pymatflow.git
```

## Usage:
There are two ways to use pymatflow:
* use Pymatflow as a python library and build your own calculation.
* use scripts provided by Pymatflow directly from command line.

The document for use Pymatflow as a python library is not ready at present! So we provide an example of using the script provided byPymatflow to show you how Pymatflow works. 
We take the scf calculation of one H2O molecule using Quantum ESPRESSO calculator as an example.
First, preparing the following structure file `h2o.xyz`(xyz format, with the second line specifying the cell parameter):
```
3
cell: 10.0   0.0   0.0  |   0.0   10.0   0.0000000  |   0.0000000   0.0000000   10.0
H    0.000000    2.599887    5.074274
H    0.000000    3.543756    3.781971
O    0.000000    2.585050    4.077973
```
make sure you put the xyz structure file and the needed pseudopotential file for H and O in the same directory and in that dir, do:
```
$ qe-scf.py -f h2o.xyz -d h2o-scf --ecutwfc 60 --conv-thr 1.0e-5 --kpoints-option gamma --runopt=genrun --auto=0
```
The above command assumes you want to run the pw.x directly on your computer, and all the calculation happens inside directory`h2o-scf`.
Usually, we do DFT calculations on powerful cluster, and pymatflow can also automatically submit the job. However, only pbs job manager and modified version of slurm on TianHe II is supported now!
For doing the calculation on the server, you should run the command on your server, and modifying the value of `--runopt` and `--auto` as follows:
```
$ qe-scf.py -f h2o.xyz -d h2o-scf --ecutwfc 60 --kpoints-option gamma --runopt=gen --auto=3 --servere=pbs
```
Which is actually the default setting for value of `--runopt` and `--auto`. This way pymatgen will automatically generate the input file and submit the job the pbs job manager. The support for other type of jobnamager is  actually not hard to implement. If you have the need, you can get in touch with me through pymatflow@163.com.
## Goals
* powerful and stable calculator supports for popular Ab initio programs
* fully featured results analysis system
* user friendly structure builder
* automated workflow for high throughput calculation[including database handling]
  * automatical choice of appropriate parameters
  * automatical correction of error during running
* statistical analysis for properties, basic machine learning module
* support for system visulization by third-party molecular visualizer like Avogadro, Vipster, paraview, vtk, vmd, xcrysden, vesta, blender
* support for third-party structure generation tools, like mbuild, atomsk, msmbuilder,
* using spgroup to get space group of crystals.

## Documents & Resources
[Official Document(ongoing!)](https://pymatflow.readthedocs.io/en/latest/)


## License
```
MIT License

Copyright (c) 2019 DeqiTang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

